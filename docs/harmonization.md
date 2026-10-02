# Harmonization layer

The reason this pipeline exists: the callers disagree on *names* even when they
agree on *biology*, and sometimes they disagree on biology too. Harmonization
separates those two cases.

## The table: `assets/harmonization.tsv`

| column | meaning |
|--------|---------|
| `gene` | gene symbol, e.g. `CYP2D6`, `GSTM1` |
| `tool` | `pypgx`, `aldy`, `cyrius`, `stellarpgx`, `pharmcat`, `optitype`, or `*` for any |
| `tool_allele` | the allele string as the tool emits it (after syntax normalization) |
| `canonical_allele` | the agreed name we map it to |
| `relation` | `synonym`, `subsumed_by`, `deletion`, `duplication`, `unsupported` |
| `note` | free text for the curator |
| `curator` | who added the row |
| `date` | when |

### `relation` values

- **synonym** — pure naming difference, same allele. Safe to collapse.
- **deletion** — a whole-gene deletion / null, however the tool spelled it.
- **duplication** — a gene duplication; the specific copy number is kept in the
  `copy_number` column of the long table, the allele is normalized to `*NxN`.
- **subsumed_by** — the tool lacks the precise allele and reports a parent/less
  specific allele (e.g. calls `*1` where a richer DB would call `*38`). Mapping
  it to the richer allele's parent keeps the comparison honest instead of
  silently "agreeing" on `*1`.
- **unsupported** — the gene/allele isn't implemented by this tool; used to
  document gaps so a missing call isn't mistaken for a reference call.

## Workflow for growing the table

1. Run the pipeline. Open each `*.unmapped.tsv` — every row is an allele spelling
   no rule or table entry covered.
2. For each, decide the `relation` and `canonical_allele` and add a row.
   Re-running re-harmonizes with no code change.
3. Review `*.discordance.tsv`. A `TIE`/`NONE` that is purely a naming artifact
   means a missing `synonym` row; a real biological disagreement stays as
   discordance and is the thing worth a human's attention.

## The GSTM1 example

Your example — one tool says `del`, another `*2`, another `-`:

```tsv
GSTM1	stellarpgx	del	*0	deletion	…
GSTM1	*	-	*0	deletion	…
GSTM1	*	*2	*0	deletion	some tools label the deletion as *2 …
GSTM1	*	present	*1	synonym	…
```

After this, all three spellings collapse to `*0` (null) and a genuine
`*0/*0` vs `*0/*1` disagreement between tools surfaces as discordance rather than
being masked by the naming difference.

## Syntax rules (code, not table)

Applied before the table lookup, in `bin/harmonize.py`:

- deletion tokens `{del, deletion, -, 0, *0, null, none}` → gene default
  (`*5` for CYP2D6, `*0` otherwise; extend `DELETION_CANONICAL`).
- duplication `*Nx<digit>` / `*N×<digit>` → `*NxN`.
- unicode asterisk `∗` → `*`.

Keep the table for anything that needs a human decision; keep the syntax rules
for things that are always mechanical.
