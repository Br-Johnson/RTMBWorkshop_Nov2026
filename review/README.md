# Review of the workshop materials

- **[`RTMB_workshop_review.md`](RTMB_workshop_review.md)** is the review: what is missing, what is wrong or poor practice, and what more modern approaches exist. Start with its section 1.
- **[`CHANGES.md`](CHANGES.md)** lists what was then changed in the workshop files (`modules/` and `examples/`), which findings each change answers, and what was left alone.
- **[`verification/`](verification/)** holds the R scripts that back every claim tagged `[run]` in the review, and **[`verification/output/`](verification/output/)** holds the console output from running them (R 4.3.3, RTMB 2.0, TMB 1.9.25). The scripts quote the slides as they were before the changes.

To rerun everything from the project root:

```sh
Rscript review/verification/run_all.R
```

The review describes the files at commit `2321223`. The corrections came afterwards, in the same pull request.
