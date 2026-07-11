# Planet Nine paper-watch runbook (for the cron agent)

You are an automated agent run on a schedule. Your ONE job: find **new** Planet
Nine / distant-perturber papers on arXiv that this repo has not yet recorded,
and write them into a findings file. Do not refactor code, do not open random
PRs, do not do anything else. Follow these steps **exactly and in order**.

Today's date is given to you in your context as `currentDate` (format
`YYYY-MM-DD`). Call it `$TODAY`. Use it verbatim where the steps say `$TODAY`.

The repo root is `/home/meawoppl/repos/planet9`. Run every command from there.

---

## Output contract (what "done" means)

By the end of a run you must have:

1. Appended every arXiv ID you evaluated this run to `paper_watch/ledger.txt`
   (so the next run never re-reports them).
2. IF and only if you found ≥1 genuinely-new relevant paper: created
   `paper_watch/findings-$TODAY.md` listing them (template in Step 6), AND
   opened one GitHub issue (Step 7).
3. Printed a one-line summary: `paper-watch $TODAY: N new, M evaluated`.

If you found **zero** new papers: do NOT create a findings file, do NOT open an
issue. Just append to the ledger and print the summary. Silence is success.

---

## Step 0 — Load the tools you need

The web tool may be deferred. Before using it, run:

```
ToolSearch  with query:  select:WebFetch
```

You will use `WebFetch` (NOT `curl` — raw network from `Bash` is sandboxed and
returns empty). You will use `Bash` for `grep`/`git`/`gh`/file writes.

---

## Step 1 — Build the "already known" set

Run this and keep the output; these are IDs you must NOT report as new:

```
cat paper_watch/ledger.txt | sort -u > /tmp/known.txt
wc -l /tmp/known.txt
```

`paper_watch/ledger.txt` already contains every arXiv ID the repo covers plus
every ID seen by past runs. Treat `/tmp/known.txt` as the source of truth for
"already known". (It was seeded from `crates/*/src/*.rs` and
`LITERATURE_GAPS.md`; you do not need to re-scan those.)

---

## Step 2 — Pull recent arXiv papers

Call `WebFetch` once per URL below. For **each** call use this exact prompt:

> "List EVERY entry in this Atom feed as one line each, newest first, in the
>  format: `YYYY-MM-DD | arXivID | first-author-surname | title`. The arXivID is
>  the part after `/abs/` without the version suffix (e.g. `2507.22297`). Do not
>  summarize, do not skip any entry, do not add commentary."

URLs (call all five):

1. `http://export.arxiv.org/api/query?search_query=all:%22Planet+Nine%22&start=0&max_results=40&sortBy=submittedDate&sortOrder=descending`
2. `http://export.arxiv.org/api/query?search_query=all:%22Planet+X%22+AND+cat:astro-ph.EP&start=0&max_results=30&sortBy=submittedDate&sortOrder=descending`
3. `http://export.arxiv.org/api/query?search_query=abs:%22extreme+trans-Neptunian%22&start=0&max_results=30&sortBy=submittedDate&sortOrder=descending`
4. `http://export.arxiv.org/api/query?search_query=abs:sednoid+OR+abs:sednoids&start=0&max_results=20&sortBy=submittedDate&sortOrder=descending`
5. `http://export.arxiv.org/api/query?search_query=abs:%22trans-Neptunian%22+AND+abs:clustering&start=0&max_results=20&sortBy=submittedDate&sortOrder=descending`

Pool all the lines from all five calls into one list. Remove exact duplicate
arXivIDs (the same paper appears in several queries — that is expected).

If a `WebFetch` call errors or returns nothing, retry it ONCE. If it still
fails, note it in the summary (`web fetch failed for query N`) and continue with
whatever you got. Never invent entries.

---

## Step 3 — Keep only RELEVANT papers

For each pooled entry, decide KEEP or DROP from its title/author line.

KEEP if it is about any of:
- Planet Nine, Planet X, Planet Y, a distant/undiscovered Solar-System planet or
  perturber, a hidden/dark Solar-System companion.
- Extreme trans-Neptunian objects (ETNOs), sednoids, detached Kuiper-belt
  objects, their orbital **clustering / alignment / precession**.
- Searches for a distant Solar-System planet (optical, far-IR/thermal, parallax,
  shift-stack, occultation, ephemeris/Cassini ranging, microlensing).
- Alternatives to Planet Nine for the clustering (self-gravitating disk, MOND /
  external field, stellar flyby, rogue/primordial planet, observational bias).

DROP (these are false positives that recur — be strict):
- **Exoplanets** around other stars (WASP / Kepler / TESS / radial-velocity
  planet detections, "nine giant planets", hot Jupiters, habitability).
- Generic Kuiper-belt / Centaur / comet papers with no distant-perturber angle.
- Interstellar objects (e.g. 3I/ATLAS) unless explicitly about Planet Nine.
- Anything clearly off-topic (machine learning, instrumentation-only, cosmology,
  galaxies) — e.g. a transformer/card-game paper is NOT relevant.

When unsure, KEEP it (a human reviews the findings; a false keep is cheap, a
missed paper is not).

---

## Step 4 — Drop the ones we already know

For each KEPT entry, in order:

1. **By ID.** If its arXivID is in `/tmp/known.txt`, DROP it (already known).
   Quick check: `grep -qxF "<arXivID>" /tmp/known.txt && echo KNOWN || echo NEW`.

2. **By author+title (catches papers implemented without the ID written down).**
   For entries still NEW after step 1, run a repo grep using the author surname
   and TWO distinctive title words, e.g.:
   `grep -rilE "russell|albedo" crates/ LITERATURE_GAPS.md`
   - If that finds a file whose content clearly matches the same paper, this is
     **AMBIGUOUS**: list it in the report under "Probably already covered —
     human check" instead of "New", and still add its ID to the ledger.
   - If nothing matches, it is genuinely **NEW**.

Keep two lists: `NEW` and `AMBIGUOUS`.

Hard rule against a known trap: **the only date that matters is the arXiv
submission date in the feed.** Ignore press coverage, anniversaries, or "news"
framing. A 2017 paper resurfacing in 2026 news is NOT new.

---

## Step 5 — Record every evaluated ID in the ledger (idempotency)

Append EVERY arXivID you evaluated this run (KEPT or DROPPED, NEW, AMBIGUOUS, or
already-known) to the ledger, then de-duplicate it:

```
# write all evaluated IDs (one per line) to /tmp/seen.txt first, then:
cat paper_watch/ledger.txt /tmp/seen.txt | sort -u -o paper_watch/ledger.txt
```

This guarantees the next run does not re-surface the same papers.

---

## Step 6 — Write the findings file (only if NEW or AMBIGUOUS is non-empty)

Create `paper_watch/findings-$TODAY.md` with exactly this shape:

```
# Planet Nine paper watch — $TODAY

## New
- **<arXivID>** (<date>) <surname> — <title>
  https://arxiv.org/abs/<arXivID>
  why relevant: <one short clause>

## Probably already covered — human check
- **<arXivID>** (<date>) <surname> — <title>  (matches: <repo file>)
  https://arxiv.org/abs/<arXivID>
```

Omit a section if its list is empty. One bullet per paper. No other prose.

---

## Step 7 — Open ONE GitHub issue (only if the `New` list is non-empty)

```
gh issue create \
  --title "Planet Nine paper watch: N new ($TODAY)" \
  --body-file paper_watch/findings-$TODAY.md
```

Replace `N` with the count in the `New` section. Do this at most once per run.
If `gh` fails, leave the findings file in place and report the failure in the
summary; do not retry more than once.

---

## Step 8 — Commit the ledger and findings

Develop on a branch, never commit to `main` (repo rule). Add files individually
(never `git add -A`). Branch name must be prefixed `meawoppl/`.

```
git checkout -b meawoppl/paper-watch-$TODAY
git add paper_watch/ledger.txt
# add the findings file too, if you created one:
git add paper_watch/findings-$TODAY.md   # skip this line if no findings file
git commit -m "Paper watch $TODAY: N new"   # title <= 10 words, NO AI attribution
git push -u origin meawoppl/paper-watch-$TODAY
```

Do NOT add any attribution to Claude/AI in the commit message or issue body
(hard repo rule). If there were zero new papers you still commit the ledger
update (the only change), with message `Paper watch $TODAY: 0 new`.

Then print: `paper-watch $TODAY: N new, M evaluated`.

---

## Do-NOT list (guardrails)

- Do NOT use `curl`/`wget` for arXiv — use `WebFetch`. Raw network is sandboxed.
- Do NOT implement crates, edit scene/figure code, or touch anything outside
  `paper_watch/`. Reproducing a paper is a human decision, not yours.
- Do NOT commit to `main`. Do NOT `git add -A`. Do NOT add AI attribution.
- Do NOT report a paper whose ID is already in the ledger.
- Do NOT invent arXiv IDs, dates, or titles. If a feed is empty, say so.
- Do NOT open more than one issue per run.

## Appendix — quick reference

- Ledger (source of truth for "known"): `paper_watch/ledger.txt`
- Findings per run: `paper_watch/findings-<date>.md`
- arXiv ID regex: `[0-9]{4}\.[0-9]{4,5}`
- The five feed URLs are in Step 2. arXiv groups months as `YYMM` in the ID:
  `2506.xxxxx` = June 2025, `2606.xxxxx` = June 2026.
