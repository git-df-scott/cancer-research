# Retiring `lymphoma-tce-architecture`

Decision: make `claude/lung-cancer-research-vl7oq3` the repository default and retire
`lymphoma-tce-architecture`, leaving one branch.

## Safety, established before anything is deleted

```
$ git merge-base --is-ancestor origin/lymphoma-tce-architecture origin/claude/lung-cancer-research-vl7oq3
$ git log --oneline origin/lymphoma-tce-architecture ^origin/claude/lung-cancer-research-vl7oq3
(no output)
```

Its tip, `30ac20b`, is an **ancestor** of the branch that will become default. Every commit it
contains is already reachable from the new default, so deleting the branch ref loses no history and
no archive tag or backup branch is required — the commits cannot be garbage-collected while the new
default exists.

Restoring it at any time, should that ever be wanted:

```bash
git branch lymphoma-tce-architecture 30ac20b
git push origin lymphoma-tce-architecture
```

## The one step I cannot perform

**Changing a repository's default branch is a GitHub settings operation.** There is no MCP tool for
it in this session and no `gh` CLI, so it needs one manual action:

> **Settings → General → Default branch → switch to `claude/lung-cancer-research-vl7oq3`**
> (https://github.com/git-df-scott/cancer-research/settings)

GitHub refuses to delete the default branch, so this must happen first. Attempting deletion before
the switch would fail rather than cause damage.

## Sequence

1. **You:** switch the default branch in Settings.
2. **Me:** delete `lymphoma-tce-architecture` from the remote.
3. **Consequence to expect:** PR #4 has `lymphoma-tce-architecture` as its base. Deleting a base
   branch closes the PR. That is correct here rather than a loss — once its head branch *is* the
   default, there is nothing left to merge. The PR's full description and history remain readable
   at https://github.com/git-df-scott/cancer-research/pull/4.

## Note on what this does not change

Retiring the branch is a repository-hygiene action. It does not merge, validate or promote any
result. The calibration remains classified **B — FIT BUT NOT VALIDATED**, and experiments R and T
remain blocked by `calibration.py` regardless of which branch is default.
