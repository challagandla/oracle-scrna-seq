# GitHub Push Instructions

## One-time setup on Ubuntu

```bash
# Configure git identity (once per machine)
git config --global user.name  "Your Name"
git config --global user.email "you@example.com"

# SSH key for GitHub (recommended)
ssh-keygen -t ed25519 -C "you@example.com"
cat ~/.ssh/id_ed25519.pub   # copy this to GitHub → Settings → SSH Keys
```

## Initialize and push this repo

```bash
# Navigate to repo root
cd scrna-pipeline

# Initialize git
git init
git branch -M main

# Add all tracked files (data/ results/ logs/ are in .gitignore)
git add .
git status            # review what will be committed

# First commit
git commit -m "Initial commit: scRNA-seq Snakemake pipeline"

# Create repo on GitHub (via CLI, requires gh tool)
gh repo create scrna-pipeline --public --description "End-to-end scRNA-seq pipeline (Snakemake + conda)"

# Or create manually at https://github.com/new, then:
git remote add origin git@github.com:<your-username>/scrna-pipeline.git

# Push
git push -u origin main
```

## After adding GitHub remote — update README badge

Edit `README.md` line 3 and replace `your-org` with your GitHub username:
```
[![CI](https://github.com/<your-username>/scrna-pipeline/actions/workflows/ci.yml/badge.svg)]
```

## Ongoing workflow

```bash
# Create a feature branch
git checkout -b feature/add-ambient-rna-step

# Make changes, then commit
git add scripts/ambient_rna.py workflow/rules/ambient.smk
git commit -m "Add: ambient RNA correction step (SoupX)"

# Push branch and open PR
git push origin feature/add-ambient-rna-step
gh pr create --title "Add ambient RNA correction" --body "Uses SoupX..."
```

## Tagging a release

```bash
git tag -a v1.0.0 -m "First stable release"
git push origin v1.0.0
```
