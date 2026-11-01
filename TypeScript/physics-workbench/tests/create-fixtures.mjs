// Writes synthetic browser-test fixtures only; contains no book material.
import { mkdirSync, writeFileSync } from "node:fs";
import { resolve } from "node:path";
const dir = resolve(process.argv[2] || "/tmp/mathphysics-workbench-fixtures");
mkdirSync(resolve(dir, "hidden"), { recursive: true });
writeFileSync(
  resolve(dir, "sample.pdf"),
  Buffer.from(
    "JVBERi0xLjQKMSAwIG9iago8PCAvVHlwZSAvQ2F0YWxvZyAvUGFnZXMgMiAwIFIgPj4KZW5kb2JqCjIgMCBvYmoKPDwgL1R5cGUgL1BhZ2VzIC9LaWRzIFszIDAgUl0gL0NvdW50IDEgPj4KZW5kb2JqCjMgMCBvYmoKPDwgL1R5cGUgL1BhZ2UgL1BhcmVudCAyIDAgUiAvTWVkaWFCb3ggWzAgMCA2MTIgNzkyXSAvUmVzb3VyY2VzIDw8IC9Gb250IDw8IC9GMSA0IDAgUiA+PiA+PiAvQ29udGVudHMgNSAwIFIgPj4KZW5kb2JqCjQgMCBvYmoKPDwgL1R5cGUgL0ZvbnQgL1N1YnR5cGUgL1R5cGUxIC9CYXNlRm9udCAvSGVsdmV0aWNhID4+CmVuZG9iago1IDAgb2JqCjw8IC9MZW5ndGggMTE4ID4+CnN0cmVhbQpCVCAvRjEgMjIgVGYgNTAgNzQwIFRkIChNYXRocGh5c2ljcyBQREYgcmVhZGluZyB0ZXN0KSBUaiAwIC0zNiBUZCAvRjEgMTQgVGYgKEEgbG9jYWwgUERGIHBhaXJlZCB3aXRoIHNhbXBsZS50ZXgpIFRqIEVUCmVuZHN0cmVhbQplbmRvYmoKeHJlZgowIDYKMDAwMDAwMDAwMCA2NTUzNSBmIAowMDAwMDAwMDA5IDAwMDAwIG4gCjAwMDAwMDAwNTggMDAwMDAgbiAKMDAwMDAwMDExNSAwMDAwMCBuIAowMDAwMDAwMjQxIDAwMDAwIG4gCjAwMDAwMDAzMTEgMDAwMDAgbiAKdHJhaWxlcgo8PCAvU2l6ZSA2IC9Sb290IDEgMCBSID4+CnN0YXJ0eHJlZgo0ODAKJSVFT0YK",
    "base64",
  ),
);
writeFileSync(
  resolve(dir, "sample.tex"),
  "\\section{Spinor notes}\nA readable identity:\n\\[ \\{\\gamma^\\mu, \\gamma^\\nu\\} = -2 g^{\\mu\\nu} I \\]\n",
);
writeFileSync(resolve(dir, "notes.md"), "# Test notes\n$E=mc^2$\n");
writeFileSync(resolve(dir, "hidden/.small.txt"), "Hidden file test\n");
console.log(dir);
