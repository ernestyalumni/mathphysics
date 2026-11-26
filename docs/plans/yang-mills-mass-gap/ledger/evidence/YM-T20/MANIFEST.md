# YM-T20 evidence manifest

Generated UTC 2026-11-26T08:03:49+00:00 by the strategist sync script.

Kept outside git (case text includes N1/N2 source text that is not ours to publish; raw responses stay private): the frozen case files, every per-attempt receipt with prompt hash and full model response, run logs. This manifest pins them.

| Private file | SHA-256 |
|---|---|
| `frozen-test.json` | `b90342fb05b0467c673bd7578018ad97c039d8fdfd668082b26cb6e3f06193fe` |
| `frozen-discovery.json` | `530c25ab54a20733826e418641f53c784e1e311f499ef148f1a9db291a7a97db` |
| `full-replay/run_full_replay.py` | `e5557fe4f89a60142257f73e95bcecc4295e0e064f0d700b128219a142332e30` |
| `full-replay/summarize_replay.py` | `9be9b7842a02a0959098f8cc4549cb9aa8ccf1453ba3f138241829134254ab85` |
| `full-replay/PREREGISTRATION.md` | `27b3ad1214488dd4db40832567a29766a276c3115ede72b66a4ecdd434c8c591` |
| `full-replay/PREREGISTRATION-GLM53.md` | `722f11df4947b593a073d9b805d33b14f714b07f0ced8f17851cc182d3ce1c60` |

The two scripts under `scripts/` differ from the private originals only in the lines that set the root paths.

| Route | Receipts | SHA-256 of sorted "path sha256" lines |
|---|---:|---|
| deepseek | 132 | `425226200fbc3c85bf2015e00c11bf7ae6e0dca01ed78a568460b28886179365` |
| flashx | 156 | `577c410132ce0bfabe1cd99425de86d86f5ee95799a2eff2fbf9202699e8d089` |
| glm53 | 87 | `f061c981821725139f0e5b6fe631f4ec2fef9fadc87f4b8e609e2169f7900e2a` |
| kimi | 86 | `dc91a690c81165d413e7a884c06da8c6931a342a9e3677c5335aeff289dd6389` |
| minimax | 81 | `aaab0f08bf3c83dbfb4f422a267ba7ceb2af04d055ee24cf931be0a211b0cf90` |
| qwen38 | 86 | `5157aa9341323acc2396dca70b0750cf445d07d26e720eea608fb2f29d5533b4` |

| Model | Runner state | Updated UTC |
|---|---|---|
| deepseek/deepseek-v4.1-flash-1017-free | complete | 2026-11-26T04:34:18.035552+00:00 |
| z-ai/glm-5.3-flashx | complete | 2026-11-26T08:01:23.054715+00:00 |
| z-ai/glm-5.3 | complete | 2026-11-26T05:56:43.873188+00:00 |
| moonshotai/kimi-k3.1-pro | complete | 2026-11-26T04:12:30.497318+00:00 |
| MiniMax-M3 | running | 2026-11-26T06:31:14.782674+00:00 |
| qwen/qwen3.8-max | complete | 2026-11-26T08:00:34.616723+00:00 |
