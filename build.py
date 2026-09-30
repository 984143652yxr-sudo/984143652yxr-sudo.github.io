"""Build a small Markdown publication into _site for GitHub Pages."""
import argparse
import html
import json
import re
import shutil
from pathlib import Path
from string import Template
import markdown

ROOT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument('--base-path', default=None, help='Empty for username.github.io; /repo for a project site')
args = parser.parse_args()
config = json.loads((ROOT / 'site.json').read_text())
base = (args.base_path if args.base_path is not None else config['base_path']).rstrip('/')
if base and not re.fullmatch(r'/[A-Za-z0-9_.-]+', base):
    raise ValueError('base_path must be empty or /repository-name')
out = ROOT / '_site'
out.mkdir(exist_ok=True)
# Rebuild only our generated output, never touch the content sources.
for child in out.iterdir():
    if child.is_dir(): shutil.rmtree(child)
    else: child.unlink()
shutil.copytree(ROOT / 'assets', out / 'assets')
template = Template((ROOT / 'templates/page.html').read_text())
esc = html.escape

def render_md(text):
    text = text.replace('{{base}}', base)
    rendered = markdown.markdown(text,
        extensions=['fenced_code', 'tables', 'toc', 'md_in_html', 'pymdownx.arithmatex'],
        extension_configs={'pymdownx.arithmatex': {'generic': True}})
    return rendered.replace('<table>', '<div class="table-wrap"><table>').replace('</table>', '</table></div>')

def page(route, title, content, active, description=None):
    dest = out / route / 'index.html'
    dest.parent.mkdir(parents=True, exist_ok=True)
    nav = ''.join(f'<a href="{base}{url}"' + (' aria-current="page"' if key == active else '') + f'>{label}</a>'
                  for key, url, label in [('home','/','Home'),('publications','/publications/','Publications'),('topics','/topics/','Study Topics'),('diary','/diary/','Diary')])
    github = config.get('github', '')
    if github and not re.fullmatch(r'https://github\.com/[A-Za-z0-9-]+/?', github):
        raise ValueError('github must be a GitHub profile URL')
    dest.write_text(template.substitute(title=esc(title), name=esc(config['name']),
        initials=esc(config['initials']), description=esc(description or config['description']),
        tagline=esc(config['tagline']), base=base, navigation=nav, content=content,
        math_script='<script defer src="https://cdn.jsdelivr.net/npm/mathjax@4.0.0/tex-mml-chtml.js"></script>' if 'class="arithmatex"' in content else '',
        affiliation=esc(config.get('affiliation', '')),
        email=f'<a href="mailto:{esc(config["email"])}">{esc(config["email"])}</a>' if config.get('email') else '',
        github=f'<a href="{esc(github)}">GitHub ↗</a>' if github else ''), encoding='utf-8')

posts = []
for path in (ROOT / 'content/posts').glob('*.md'):
    # Metadata is a JSON object between the first two --- lines.
    _, metadata, body = path.read_text().split('---', 2)
    p = json.loads(metadata)
    if not re.fullmatch(r'[a-z0-9-]+', p['slug']): raise ValueError('Invalid post slug')
    p['body'] = body
    posts.append(p)
posts.sort(key=lambda p: p['date'], reverse=True)
if len({p['slug'] for p in posts}) != len(posts): raise ValueError('Duplicate post slug')

def entry(p, featured=False):
    return f'''<article class="{'featured' if featured else 'entry'}">
    <div class="meta"><time datetime="{esc(p['date'])}">{esc(p['date_label'])}</time><span class="badge">{esc(p['status'])}</span></div>
    <h2><a href="{base}/diary/{p['slug']}/">{esc(p['title'])}</a></h2><p>{esc(p['summary'])}</p>
    <a class="text-link" href="{base}/diary/{p['slug']}/">Read the study note <span aria-hidden="true">→</span></a></article>'''

home = render_md((ROOT / 'content/home.md').read_text())
# Keep one publication source for both the home section and the full page.
publications_source = (ROOT / 'content/publications.md').read_text()
publication_body = publications_source.split('\n', 1)[1].strip()
home = home.replace('<!-- publications -->', render_md(publication_body))
page('', 'Home', home, 'home')
page('publications', 'Publications', render_md(publications_source), 'publications')
topics_content = render_md((ROOT / 'content/topics.md').read_text())
page('topics', 'Study Topics', topics_content, 'topics')
page('research', 'Study Topics', topics_content, 'topics')
page('diary', 'Diary', '<h1>Diary</h1><p class="lead">Dated notes on papers, code, and ongoing analyses.</p>' + ''.join(entry(p) for p in posts), 'diary')
for p in posts:
    header = f'''<a class="back" href="{base}/diary/">← All study notes</a><header class="article-header">
    <h1>{esc(p['title'])}</h1><p class="lead">{esc(p['summary'])}</p>
    <div class="meta"><time datetime="{esc(p['date'])}">{esc(p['date_label'])}</time><span class="badge">{esc(p['status'])}</span><span>{esc(p['reading_time'])}</span></div></header>'''
    page('diary/' + p['slug'], p['title'], header + '<article class="prose">' + render_md(p['body']) + '</article>', 'diary', p['summary'])
page('404-page', 'Page not found', f'<h1>Page not found</h1><p>This note may have moved.</p><a href="{base}/diary/">Browse the study diary →</a>', '')
shutil.copyfile(out / '404-page/index.html', out / '404.html')
shutil.rmtree(out / '404-page')
(out / '.nojekyll').touch()
print(f'Built {len(posts)} study notes and 4 main pages in {out}')
