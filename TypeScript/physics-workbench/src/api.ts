import type { Entry, Root, Kind } from "./model";
export interface Listing {
  root: string;
  path: string;
  entries: Entry[];
  truncated: boolean;
  visited: number;
}
export interface TextFile {
  root: string;
  path: string;
  text: string;
  kind: Kind;
  paired: { path: string; kind: Kind } | null;
  modified_secs: number;
}
export async function request<T>(
  path: string,
  options: RequestInit = {},
): Promise<T> {
  const response = await fetch(path, options);
  const data = await response.json();
  if (!response.ok)
    throw new Error(
      typeof data.error === "string"
        ? data.error
        : typeof data.error?.message === "string"
          ? data.error.message
          : `Request failed (${response.status})`,
    );
  return data as T;
}
export const query = (root: string, path: string) =>
  new URLSearchParams({ root, path });
export const getRoots = () => request<{ roots: Root[] }>("/api/roots");
export const getFiles = (
  root: string,
  path: string,
  q: string,
  signal?: AbortSignal,
) =>
  request<Listing>("/api/files?" + new URLSearchParams({ root, path, q }), {
    signal,
  });
export const getText = (root: string, path: string) =>
  request<TextFile>("/api/text?" + query(root, path));
export const pdfUrl = (root: string, path: string) =>
  "/api/pdf?" + query(root, path);
export const runTrace = (labels: string[]) =>
  request<Record<string, unknown>>("/api/tool", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      command: "gamma_trace",
      labels,
      domain: "formal_symmetric_h",
      budget: 10395,
    }),
  });
