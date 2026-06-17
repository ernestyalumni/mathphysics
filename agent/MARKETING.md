# MARKETING.md — Short-Form Video & Social Media Content

Last updated: 2026-06-16

This file is TMPAI's content calendar, draft index, and publication ledger. All drafting is allowed; **publishing requires Ernest's explicit approval.**

## Why TMPAI does marketing

The mission isn't just "do physics in private." It's to **build an audience**, **recruit collaborators**, **earn attention for verified results**, and over time make the discovery engine itself a public artifact. Short-form video and social media are the cheapest, highest-leverage channels for that.

Without communication, the engine is a hobby. With it, it can fund itself and attract serious collaborators.

## Truth rules (non-negotiable)

1. **No clickbait that misstates the result.** "I solved gauge theory" when you derived Parke-Taylor is not allowed.
2. **Re-derivations are labelled.** If Srednicki ch. 48 has it, the script says "Here's how Srednicki proves..." not "Here's a new result."
3. **New results require verifier output.** Description / show-notes link to `discoveries/YYYY-MM-DD-<slug>.md` and the Cadabra2 script.
4. **No fabricated citations.** Every paper named is one Ernest has read.
5. **No politics, no rage-bait, no influencer drama.** Stay in lane: math, physics, AI tools for them, occasional dev-process content.

## Channels

| Channel | Format | Length | Status |
|---------|--------|--------|--------|
| YouTube Shorts | Vertical, 60-90s | 60-90s | Planned |
| TikTok | Vertical, 60-90s | 60-90s | Planned (mirror of Shorts) |
| Instagram Reels | Vertical, 60-90s | 60-90s | Planned (mirror of Shorts) |
| YouTube long-form | 16:9, narrative | 8-20m | Planned |
| X / Twitter | Threads + clips | 4-10 tweets | Planned |
| Blog / personal site | Long-form text | 2k-6k words | Planned |
| arXiv | Papers | full | Already in pipeline (MHV paper v2) |

## Style references

- **3Blue1Brown** — Manim, narrative + clear math, slow zoom on the central equation.
- **Welch Labs** — visual intuition + honest engagement with the actual math, willing to show the gnarly bits.
- **Visualizing AI / Two Minute Papers** — short, dense, one idea per video, results-led.
- **Numberphile** — interview-style, personality-led; less reproducible but works for one-off collabs.
- **DeepMind / OpenAI research videos** — for the engine pieces; very visual, big-frame narration.

## Content pipeline

```
Verified result (from a loop or session) or shippable re-derivation
   ↓
TMPAI drafts:
  - 60-90s short script    → marketing/drafts/YYYY-MM-DD-<slug>-short.md
  - 5-7m long outline      → marketing/drafts/YYYY-MM-DD-<slug>-long.md
  - 4-7 tweet thread       → marketing/drafts/YYYY-MM-DD-<slug>-thread.md
  - 1 blog post outline    → marketing/drafts/YYYY-MM-DD-<slug>-blog.md
   ↓
Status: DRAFT — awaiting Ernest's approval
   ↓
Ernest reviews. If approved:
  - Ernest publishes (or TMPAI publishes when explicit posting is enabled per channel)
  - Move draft to marketing/published/YYYY-MM-DD-<slug>.md
  - Update the published index in this file
```

## Draft file template (`marketing/drafts/YYYY-MM-DD-<slug>-short.md`)

```markdown
# <slug> — 60-90s short

**Status:** DRAFT — awaiting Ernest's approval
**Channel:** YouTube Shorts / TikTok / Reels (vertical 9:16)
**Source result:** <link to discoveries/YYYY-MM-DD-<slug>.md OR cite the paper section>
**Hook archetype:** <e.g. "weird symmetry", "computer-verified", "the line that took 50 years">

## 0-3s — Hook
<single-sentence visual + spoken hook>

## 3-15s — Setup
<establish the problem in plain English, one screen of context>

## 15-60s — The math beat
<the actual result, with the central equation visible the whole time>

## 60-80s — The "why care" beat
<one sentence connecting to a bigger thread — gauge invariance, gravity, amplitude program>

## 80-90s — Closing hook
<one sentence + CTA: "full derivation on the channel" or "Cadabra2 script in the repo">

## Visual notes
- <Manim shot list, or whiteboard cues, or screencap timestamps>

## Music / pacing
- <reference track, BPM band, where to cut>

## Truth check
- Re-derivation? Cite source: <Srednicki ch. N | arXiv:NNNN.NNNNN>
- Verifier output: <link to Cadabra2 script + residual>
- Anything ambiguous I'm rounding? <answer or say none>
```

## Content calendar (target cadence once Stage 1 hits)

- **Weekly short** — 1 verified-result or re-derivation short per week.
- **Monthly long-form** — 1 deep dive per month.
- **Weekly thread** — 1 results thread per week (often tied to the short).
- **Quarterly blog** — 1 long-form blog post per quarter, ties everything together.

Don't enforce cadence before the first 3 shorts ship with Ernest's approval.

## Draft index

(Empty at seeding. TMPAI updates this whenever a new draft is created.)

| Date | Slug | Channel | Status | Source result |
|------|------|---------|--------|---------------|
| – | – | – | – | – |

## Published index

(Empty at seeding. TMPAI updates this once Ernest publishes.)

| Date published | Channel | Slug | URL | Notes |
|----------------|---------|------|-----|-------|
| – | – | – | – | – |

## Anti-pattern catalog

- **"AI just discovered ____."** No. AI re-derived a thing, or AI proposed and a verifier checked. Be specific.
- **Pretending Manim animations are simulations.** They're animations. Say so if asked.
- **Citing a paper the team hasn't read.** Don't.
- **Cross-posting before Ernest approved each channel.** Approval is per-channel, not blanket.
- **Trend-chasing.** Drop the trend if it would require misstating the math.

## When to ship vs hold

Ship a short when:

- It's a verified result + a clean visual beat exists + the truth check is green.
- It's a clearly-labelled re-derivation of a beautiful classical result.
- It's a process / behind-the-scenes piece about the engine itself (clearly framed as process).

Hold a short when:

- The verifier was inconclusive.
- The novelty claim hasn't been literature-checked.
- Ernest hasn't approved yet.
- The hook would only work by misstating the math.

## Cross-references

- `../TMP-MISSION.md` — pillar 4.
- `DISCOVERY-LOOPS.md` — content sources.
- `RESEARCH-AGENDA.md` — current threads and their shorts potential.
