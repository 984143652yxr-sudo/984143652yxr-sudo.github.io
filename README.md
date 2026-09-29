# Catherine Xinrui Yu — academic homepage

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
