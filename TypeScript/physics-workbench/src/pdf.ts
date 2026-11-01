import {
  getDocument,
  GlobalWorkerOptions,
  TextLayer,
} from "pdfjs-dist/legacy/build/pdf.mjs";
import type {
  PDFDocumentProxy,
  RenderTask,
} from "pdfjs-dist/types/src/display/api";
import workerUrl from "pdfjs-dist/legacy/build/pdf.worker.min.mjs?url";
GlobalWorkerOptions.workerSrc = workerUrl;
export interface PdfView {
  goTo(page: number): void;
  setZoom(zoom: number): void;
  dispose(): void;
}
/** One live page at a time, backed by local PDF.js assets and a selectable text layer. */
export function mountPdf(
  host: HTMLElement,
  url: string,
  initialPage: number,
  onPage: (page: number, total: number) => void,
): PdfView {
  host.className = "pdf-view";
  host.setAttribute("aria-label", "PDF page viewer");
  const status = document.createElement("p");
  status.className = "pdf-loading";
  status.setAttribute("role", "status");
  status.textContent = "Opening PDF…";
  host.append(status);
  const task = getDocument({
    url,
    cMapUrl: "/pdfjs/cmaps/",
    cMapPacked: true,
    standardFontDataUrl: "/pdfjs/standard_fonts/",
    wasmUrl: "/pdfjs/wasm/",
    enableXfa: false,
    disableAutoFetch: true,
  });
  let doc: PDFDocumentProxy | undefined;
  let disposed = false;
  let generation = 0;
  let pageNumber = Math.max(1, initialPage || 1);
  let zoom = 1;
  let renderTask: RenderTask | undefined;
  let textLayer: TextLayer | undefined;
  async function draw(requested: number) {
    const token = ++generation;
    if (!doc || disposed) return;
    pageNumber = Math.max(
      1,
      Math.min(doc.numPages, Math.floor(requested) || 1),
    );
    renderTask?.cancel();
    textLayer?.cancel();
    status.textContent = "Rendering page " + pageNumber + "…";
    host.replaceChildren(status);
    try {
      const page = await doc.getPage(pageNumber);
      if (disposed || token !== generation) return;
      const natural = page.getViewport({ scale: 1 });
      if (
        !Number.isFinite(natural.width) ||
        !Number.isFinite(natural.height) ||
        natural.width <= 0 ||
        natural.height <= 0
      )
        throw new Error("This PDF has invalid page dimensions.");
      const ratio = Math.min(window.devicePixelRatio || 1, 2);
      let scale = Math.max(
        0.2,
        Math.min(
          3,
          (Math.max(260, host.clientWidth - 40) / natural.width) * zoom,
        ),
      );
      const pixels =
        natural.width * natural.height * scale * scale * ratio * ratio;
      if (pixels > 16_000_000) scale *= Math.sqrt(16_000_000 / pixels);
      const viewport = page.getViewport({ scale });
      const paper = document.createElement("div");
      paper.className = "pdf-page";
      paper.style.width = viewport.width + "px";
      paper.style.height = viewport.height + "px";
      paper.style.setProperty(
        "--total-scale-factor",
        String(scale * page.userUnit),
      );
      paper.dataset.page = String(pageNumber);
      const canvas = document.createElement("canvas");
      canvas.width = Math.floor(viewport.width * ratio);
      canvas.height = Math.floor(viewport.height * ratio);
      canvas.style.width = viewport.width + "px";
      canvas.style.height = viewport.height + "px";
      canvas.setAttribute("aria-label", "Rendered PDF page " + pageNumber);
      const text = document.createElement("div");
      text.className = "textLayer";
      paper.append(canvas, text);
      host.replaceChildren(paper);
      renderTask = page.render({
        canvas,
        viewport,
        transform: ratio !== 1 ? [ratio, 0, 0, ratio, 0, 0] : undefined,
        background: "white",
      });
      await renderTask.promise;
      if (disposed || token !== generation) return;
      textLayer = new TextLayer({
        container: text,
        viewport,
        textContentSource: page.streamTextContent(),
      });
      await textLayer.render();
      if (disposed || token !== generation) return;
      paper.dataset.rendered = "true";
      onPage(pageNumber, doc.numPages);
    } catch (e) {
      if (
        disposed ||
        token !== generation ||
        (e as Error).name === "RenderingCancelledException"
      )
        return;
      status.textContent = "Could not render PDF: " + (e as Error).message;
      host.replaceChildren(status);
    }
  }
  void task.promise
    .then((pdf) => {
      if (disposed) {
        void task.destroy();
        return;
      }
      doc = pdf;
      void draw(pageNumber);
    })
    .catch((e) => {
      if (!disposed) {
        status.textContent = "Could not open PDF: " + (e as Error).message;
        host.replaceChildren(status);
      }
    });
  return {
    goTo(page) {
      pageNumber = Number.isFinite(page) ? Math.max(1, Math.floor(page)) : 1;
      void draw(pageNumber);
    },
    setZoom(value) {
      zoom = Number.isFinite(value) ? Math.min(2, Math.max(0.5, value)) : 1;
      void draw(pageNumber);
    },
    dispose() {
      disposed = true;
      generation++;
      renderTask?.cancel();
      textLayer?.cancel();
      void task.destroy().catch(() => {});
    },
  };
}
