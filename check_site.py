"""Validate generated page titles, internal links, assets, and fragment targets."""
import argparse
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

parser = argparse.ArgumentParser()
parser.add_argument('--base-path', default='')
base = parser.parse_args().base_path.rstrip('/')
root = Path(__file__).resolve().parent / '_site'

class Page(HTMLParser):
    def __init__(self, text):
        super().__init__()
        self.links, self.ids, self.h1, self.title = [], set(), 0, 0
        self.feed(text)
    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == 'h1': self.h1 += 1
        if tag == 'title': self.title += 1
        if 'id' in attrs: self.ids.add(attrs['id'])
        for key in ('href', 'src'):
            if key in attrs: self.links.append(attrs[key])

pages = {p: Page(p.read_text()) for p in root.rglob('*.html')}
errors = []
for path, page in pages.items():
    if page.h1 != 1 or page.title != 1: errors.append(f'{path}: needs exactly one h1 and title')
    for link in page.links:
        url = urlsplit(link)
        if url.scheme or url.netloc: continue
        target = unquote(url.path)
        if target.startswith('/'):
            if base and not (target == base or target.startswith(base + '/')):
                errors.append(f'{path}: missing base path in {link}')
                continue
            dest = root / target[len(base):].lstrip('/')
        elif target:
            dest = path.parent / target
        else:
            dest = path
        if dest.is_dir(): dest = dest / 'index.html'
        dest = dest.resolve()
        if not dest.exists(): errors.append(f'{path}: missing target {link}')
        elif url.fragment and dest in pages and unquote(url.fragment) not in pages[dest].ids:
            errors.append(f'{path}: missing fragment {link}')
if errors: raise SystemExit('\n'.join(errors))
print(f'Checked {len(pages)} HTML pages: internal links, fragments, assets, and headings pass.')
