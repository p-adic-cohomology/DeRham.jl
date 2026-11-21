# VarByVar performance recheck

## Setup

Compared commit A (`2f899bd`) with commit B (`034b418`) using the E2E case
`gpu_p7_varbyvar__dim_3_deg_3_fermat_2` on CPU with `use_gpu=false`,
`fastevaluation=true`, `S=[4]`, and p=7. Each commit ran in its own Julia
1.13.2 process with one Julia thread. After warm-up, measurements alternated
A/B/A/B/A/B in the same environment. Every run matched the stored expected
L-polynomial.

## Work counts

Each reduction chunk invokes the evaluation kernel once, so the chunk count is
also the evaluation count.

| Commit | Chunks / evaluations | Reduction steps | Duplicate merges |
| --- | ---: | ---: | ---: |
| A (`2f899bd`) | 13,922 | 97,835 | 4,580 |
| B (`034b418`) | 10,302 | 94,215 | 4,575 |

Chains per pass were:

- A: `[126, 70, 35, 15, 5]` for five passes, then `[792, 330, 120, 36, 8]` for five passes (50 passes total).
- B: `[126, 70, 70, 35, 35, 15, 15, 5, 5]` once; `[126, 70, 70, 35, 35, 15, 15, 4, 4]` four times; `[792, 330, 330, 120, 120, 36, 8, 8]` twice; `[792, 330, 330, 120, 36, 8, 8]` once; `[792, 330, 330, 120, 36, 36, 8, 8]` once; and `[792, 330, 330, 120, 36, 36, 7, 7]` once (85 passes total).

## Timings and disposition

| Run | A (seconds) | B (seconds) |
| ---: | ---: | ---: |
| 1 | 136.211 | 132.316 |
| 2 | 158.240 | 126.419 |
| 3 | 139.620 | 119.971 |
| Median | 139.620 | 126.419 |

B's median is 0.905 times A's median, about 9.5% faster, within the 1.15x
acceptance threshold. B also uses fewer chunks, steps, and duplicate merges.
An earlier sequential warm comparison (62.928 seconds for A and 95.343 for B)
was not reproduced by the alternating runs. The observed data does not support
a VarByVar policy change or a profiling follow-up.
