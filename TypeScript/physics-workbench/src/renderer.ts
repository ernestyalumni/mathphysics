import katex from "katex";
import { marked } from "marked";
import DOMPurify from "dompurify";

const LIMIT = 250_000;
const MACROS: Record<string, string> = {
  "\\ket": "\\left|#1\\right\\rangle",
  "\\bra": "\\left\\langle#1\\right|",
  "\\braket": "\\left\\langle#1\\right\\rangle",
  "\\mel": "\\left\\langle#1\\right\\rangle",
  "\\ev": "\\left\\langle#1\\right\\rangle",
  "\\expval": "\\left\\langle#1\\right\\rangle",
  "\\Tr": "\\operatorname{Tr}",
  "\\tr": "\\operatorname{tr}",
  "\\rank": "\\operatorname{rank}",
  "\\diag": "\\operatorname{diag}",
  "\\ii": "i",
  "\\ee": "e",
  "\\dd": "\\mathrm{d}",
  "\\dv": "\\frac{\\mathrm{d}#1}{\\mathrm{d}#2}",
  "\\pdv": "\\frac{\\partial #1}{\\partial #2}",
  "\\abs": "\\left|#1\\right|",
  "\\norm": "\\left\\|#1\\right\\|",
  "\\comm": "\\left[#1,\\,#2\\right]",
  "\\acomm": "\\left\\{#1,\\,#2\\right\\}",
  "\\order": "\\mathcal{O}\\!\\left(#1\\right)",
  "\\vb": "\\mathbf{#1}",
  "\\vu": "\\hat{\\mathbf{#1}}",
  "\\grad": "\\nabla",
  "\\curl": "\\nabla\\times",
  "\\divergence": "\\nabla\\cdot",
  "\\laplacian": "\\nabla^2",
  "\\cross": "\\times",
  "\\T": "^{\\mathsf{T}}",
  "\\R": "\\mathbb{R}",
  "\\C": "\\mathbb{C}",
  "\\N": "\\mathbb{N}",
  "\\Z": "\\mathbb{Z}",
  "\\Q": "\\mathbb{Q}",
  "\\SO": "\\mathrm{SO}(#1)",
  "\\SU": "\\mathrm{SU}(#1)",
  "\\SE": "\\mathrm{SE}(#1)",
};
const escape = (s: string): string =>
  s.replace(
    /[&<>"']/g,
    (c) =>
      ({
        "&": "&amp;",
        "<": "&lt;",
        ">": "&gt;",
        '"': "&quot;",
        "'": "&#39;",
      })[c]!,
  );
const warning = (s: string): string =>
  `<p class="renderer-warning" role="status">${escape(s)}</p>`;
const tooLarge = (s: string): string | undefined =>
  s.length > LIMIT
    ? warning(
        `Preview not rendered: ${s.length} characters exceeds the ${LIMIT}-character limit. Open the original source in an editor.`,
      )
    : undefined;

/** Render one math body (without delimiters). No persistent user macro state. */
export function renderMath(source: string, display = false): string {
  const oversized = tooLarge(source);
  if (oversized) return oversized;
  try {
    return katex.renderToString(source, {
      displayMode: display,
      trust: false,
      maxExpand: 1000,
      maxSize: 20,
      throwOnError: true,
      strict: "ignore",
      macros: { ...MACROS },
      globalGroup: false,
    });
  } catch (error) {
    const message = error instanceof Error ? error.message : String(error);
    return `<span class="render-error" role="status">${escape(message)}: <code>${escape(source)}</code></span>`;
  }
}

type Slot = { source: string; display: boolean; raw: string };
const escapedAt = (s: string, i: number): boolean => {
  let n = 0;
  while (i > 0 && s[--i] === "\\") n++;
  return n % 2 === 1;
};

/** Scan code first, before Markdown can consume TeX backslashes/underscores. */
function protect(
  source: string,
  tex: boolean,
): { source: string; slots: Slot[]; prefix: string } {
  let prefix = "MPMATHSLOT";
  while (source.includes(prefix)) prefix += "X";
  const slots: Slot[] = [];
  const slot = (body: string, display: boolean, raw: string): string => {
    slots.push({ source: body, display, raw });
    return `${prefix}${slots.length - 1}END`;
  };
  const parts: string[] = [];
  let i = 0;
  while (i < source.length) {
    const rest = source.slice(i);
    const lineStart = i === 0 || source[i - 1] === "\n";
    if (!tex && lineStart) {
      const fence = /^( {0,3})(`{3,}|~{3,})[^\n]*(?:\n|$)/.exec(rest);
      if (fence) {
        const end = new RegExp(
          `^ {0,3}${fence[2][0]}{${fence[2].length},}[ \\t]*$`,
          "gm",
        );
        end.lastIndex = fence[0].length;
        const close = end.exec(rest);
        const length = close ? close.index + close[0].length : rest.length;
        parts.push(rest.slice(0, length));
        i += length;
        continue;
      }
      if (/^(?: {4}|\t)/.test(rest)) {
        const end = rest.indexOf("\n");
        const length = end < 0 ? rest.length : end + 1;
        parts.push(rest.slice(0, length));
        i += length;
        continue;
      }
    }
    if (!tex && source[i] === "`" && !escapedAt(source, i)) {
      const ticks = /^`+/.exec(rest)![0];
      let end = source.indexOf(ticks, i + ticks.length);
      while (
        end >= 0 &&
        (source[end - 1] === "`" || source[end + ticks.length] === "`")
      )
        end = source.indexOf(ticks, end + ticks.length);
      if (end >= 0) {
        parts.push(source.slice(i, end + ticks.length));
        i = end + ticks.length;
        continue;
      }
    }
    if (!escapedAt(source, i)) {
      const env = tex ? /^\\begin\{(align\*?|equation\*?)\}/.exec(rest) : null;
      if (env) {
        const close = `\\end{${env[1]}}`;
        const end = source.indexOf(close, i + env[0].length);
        if (end >= 0) {
          const body = source.slice(i + env[0].length, end);
          parts.push(
            slot(
              env[1].startsWith("align")
                ? `\\begin{aligned}${body}\\end{aligned}`
                : body,
              true,
              source.slice(i, end + close.length),
            ),
          );
          i = end + close.length;
          continue;
        }
      }
      let open = "",
        close = "",
        display = false;
      if (rest.startsWith("$$")) {
        open = close = "$$";
        display = true;
      } else if (rest.startsWith("\\[")) {
        open = "\\[";
        close = "\\]";
        display = true;
      } else if (rest.startsWith("\\(")) {
        open = "\\(";
        close = "\\)";
      } else if (rest[0] === "$" && rest[1] && !/\s/.test(rest[1])) {
        open = close = "$";
      }
      if (open) {
        let end = source.indexOf(close, i + open.length);
        while (end >= 0 && escapedAt(source, end))
          end = source.indexOf(close, end + close.length);
        if (end >= 0) {
          const body = source.slice(i + open.length, end);
          if (
            open !== "$" ||
            (body.length > 0 &&
              !/[\s\\]$/.test(body) &&
              !body.includes("\n") &&
              !/\d/.test(source[end + 1] || ""))
          ) {
            parts.push(
              slot(body, display, source.slice(i, end + close.length)),
            );
            i = end + close.length;
            continue;
          }
        }
      }
    }
    parts.push(source[i++]);
  }
  return { source: parts.join(""), slots, prefix };
}

/** Browser DOM required. Safe returned HTML may be assigned to a preview innerHTML. */
export function renderDocument(
  source: string,
  kind: "tex" | "markdown" | "text",
): string {
  const oversized = tooLarge(source);
  if (oversized) return oversized;
  if (kind === "text")
    return `<pre class="source-text">${escape(source)}</pre>`;
  try {
    const protectedSource = protect(source, kind === "tex");
    let html: string;
    if (kind === "tex") {
      // Only standalone simple headings are interpreted; all other source stays visible.
      html = protectedSource.source
        .split("\n")
        .map((line) => {
          const heading =
            /^\s*\\(section|subsection|subsubsection|chapter|title)\*?\{([^{}]*)\}\s*$/.exec(
              line,
            );
          if (!heading) return escape(line);
          const level = {
            chapter: 1,
            title: 1,
            section: 2,
            subsection: 3,
            subsubsection: 4,
          }[heading[1]]!;
          return `<h${level}>${escape(heading[2])}</h${level}>`;
        })
        .join("\n");
      html = `<div class="tex-preview">${html}</div>`;
    } else {
      // Separate instance options: never mutate the host's global marked configuration.
      html = marked.parse(protectedSource.source, {
        async: false,
        gfm: true,
        breaks: true,
      }) as string;
    }
    // User HTML never reaches the preview unsanitized. KaTeX output is inserted later.
    const safe = DOMPurify.sanitize(html, {
      USE_PROFILES: { html: true },
      FORBID_TAGS: [
        "style",
        "script",
        "iframe",
        "object",
        "embed",
        "form",
        "input",
      ],
      FORBID_ATTR: ["style", "srcdoc"],
      ALLOW_DATA_ATTR: false,
    });
    const template = document.createElement("template");
    template.innerHTML = safe;
    // Explicit URL allowlist also excludes data:, file:, blob: and protocol-relative URLs.
    template.content
      .querySelectorAll("[href], [src], [srcset], [poster], [action]")
      .forEach((el) => {
        for (const attr of ["href", "src", "poster", "action"]) {
          const value = el.getAttribute(attr);
          if (
            value !== null &&
            !/^(?:https?:|mailto:|#|\/(?!\/)|\.\.?\/)/i.test(value.trim())
          )
            el.removeAttribute(attr);
        }
        el.removeAttribute("srcset");
      });
    const walker = document.createTreeWalker(
      template.content,
      4 /* SHOW_TEXT */,
    );
    const texts: Text[] = [];
    while (walker.nextNode()) texts.push(walker.currentNode as Text);
    const pattern = new RegExp(`${protectedSource.prefix}(\\d+)END`, "g");
    for (const node of texts) {
      const inCode = !!node.parentElement?.closest("code, pre");
      const content = node.data;
      pattern.lastIndex = 0;
      let match: RegExpExecArray | null;
      let cursor = 0;
      const fragment = document.createDocumentFragment();
      while ((match = pattern.exec(content))) {
        fragment.append(
          document.createTextNode(content.slice(cursor, match.index)),
        );
        const entry = protectedSource.slots[Number(match[1])];
        if (entry && inCode)
          fragment.append(document.createTextNode(entry.raw));
        else if (entry) {
          const math = document.createElement("template");
          math.innerHTML = renderMath(entry.source, entry.display);
          fragment.append(math.content);
        } else fragment.append(document.createTextNode(match[0]));
        cursor = match.index + match[0].length;
      }
      if (cursor) {
        fragment.append(document.createTextNode(content.slice(cursor)));
        node.replaceWith(fragment);
      }
    }
    return (
      (kind === "tex"
        ? warning(
            "Source + math preview only: complex preambles, custom document macros and full TeX layout are unsupported. Unrecognized commands remain source text.",
          )
        : "") + template.innerHTML
    );
  } catch (error) {
    return (
      warning(
        `Preview unavailable: ${error instanceof Error ? error.message : String(error)}`,
      ) + `<pre>${escape(source)}</pre>`
    );
  }
}
