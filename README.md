# Catherine Xinrui Yu — academic homepage

A plain white academic website with compact top navigation, inspired by [Lijun Wang’s homepage](https://hohoweiya.xyz/) and [Yitao Xu’s homepage](https://yitaoxu.github.io/). Original layout code and styling; no copied biography, photographs, or publications.

Includes Home (personal information and publications), Publications, Study Topics, Diary, and an initial CIGMA–OneK1K replication entry. Posts are Markdown, compiled into static HTML with Python. No browser JavaScript or backend is required.

## Current status

Configured for GitHub Pages at https://984143652yxr-sudo.github.io, with source in `984143652yxr-sudo/984143652yxr-sudo.github.io`. The academic biography is based on the author’s CV.

The first entry summarizes the local CIGMA tutorial notebooks and September 28 progress report. Its numeric checkpoint is attributed to that report, not independently recomputed. No original notebooks, data files, private filesystem paths, or personal phone number are included in the website. The academic email from the CV is displayed for professional contact.

## Preview locally

From this folder, using Python 3.10 or newer:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python build.py
.venv/bin/python check_site.py
.venv/bin/python -m http.server 8765 --bind 127.0.0.1 --directory _site
```

Open http://127.0.0.1:8765. After changing Markdown, run `build.py` again and refresh the browser. `_site/` contains generated output and is ignored by Git.

## Edit directly on GitHub

Open a source file, choose the pencil icon, make your changes, and commit to `main`. GitHub Actions rebuilds and publishes the site automatically.

- `content/home.md`: About me and research interests
- `content/publications.md`: publications and manuscript, on both Home and Publications
- `content/topics.md`: Study Topics
- `content/posts/`: diary entries
- `site.json`: name, affiliation, academic email, and GitHub profile

## Personalize

Edit `site.json` for the name, short description, and a confirmed GitHub profile URL. Edit `content/home.md` for the biography and `content/topics.md` for Study Topics. Edit `content/publications.md` once to update both the homepage publication section and the Publications page. The Biometrics citation was checked against the published article. The unpublished work is listed separately as a manuscript, using the title and author list in the CV. The full CV and phone number are not included. Colors, typography, and responsive rules live in `assets/style.css`.

For a personal site such as `username.github.io`, leave `base_path` empty. For a project site, use `/repository-name`. The publishing workflow reads this automatically from GitHub Pages metadata.

## Add a study note

1. Copy `POST-TEMPLATE.md` into `content/posts/YYYY-MM-DD-short-title.md`.
2. Edit the JSON metadata between the `---` lines. Use a unique lowercase `slug`; it becomes the URL.
3. Write the note below the metadata using Markdown.
4. Run the build and check commands above, then commit and push when ready.

Diary lists every post in reverse chronological order; the homepage remains focused on academic information and publications. The old `/research/` URL remains available as an alias of Study Topics. A `Working draft` status is a visible label, **not** a privacy control: every file in `content/posts` is rendered. Keep private drafts outside that directory.

For figures, put only shareable images in `assets/` and use `![Descriptive caption]({{base}}/assets/figure.png)`. A standalone `[TOC]` inserts a table of contents.

## Publish on GitHub Pages

Use this folder as the **root of a separate website repository**, not the parent `life general` repository. Copy its files including `.github/`, but excluding `.venv/` and `_site/`.

1. Use the confirmed account `984143652yxr-sudo` and a repository named `984143652yxr-sudo.github.io`.
2. Put this website source on the repository’s `main` branch.
3. In **Settings → Pages → Build and deployment**, select **GitHub Actions**.
4. Run **Publish study diary** from Actions, or push a change to `main`.
5. Check the deployment and open the resulting Pages URL.

The included workflow installs the pinned Markdown dependency, builds the website, checks internal links and fragments, and deploys only `_site/`.

Official references: [Create a Pages site](https://docs.github.com/en/pages/getting-started-with-github-pages/creating-a-github-pages-site) and [use custom workflows](https://docs.github.com/en/pages/getting-started-with-github-pages/using-custom-workflows-with-github-pages).
