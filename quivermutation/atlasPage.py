"""The atlas's top shapes as a page to look at, not keys to read.

`classpage` draws lines and quipus in the browser from their names; the shapes
here have no names, so they are drawn in Python from the quivers themselves --
a spring layout with a fixed seed, parallel arrows bowed apart, and each
relation written out underneath.  The page is one static file with everything
inline and follows the light and dark schemes of whatever opens it.
"""

import html
import json
import math

import networkx as nx
import polars as pl

from . import shapeAtlas
from . import shapeKeys


_STYLE = """
:root { --bg: #fbfaf7; --fg: #1d1d1b; --muted: #6b6a66; --line: #3b3a36;
        --node: #ffffff; --accent: #b4532a; --card: #f2f0ea; }
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) { --bg: #161615; --fg: #ecebe6; --muted: #9c9a93;
        --line: #cfcdc5; --node: #262624; --accent: #e08a5f; --card: #1f1f1d; }
}
:root[data-theme="dark"] { --bg: #161615; --fg: #ecebe6; --muted: #9c9a93;
        --line: #cfcdc5; --node: #262624; --accent: #e08a5f; --card: #1f1f1d; }
body { background: var(--bg); color: var(--fg); margin: 0 auto; max-width: 1100px;
       padding: 24px 16px; font: 15px/1.5 system-ui, sans-serif; }
h1 { font-size: 1.5rem; } h2 { font-size: 1.15rem; margin-top: 2rem; }
.grid { display: grid; gap: 16px; grid-template-columns: repeat(auto-fill, minmax(260px, 1fr)); }
.card { background: var(--card); border-radius: 8px; padding: 12px; overflow-wrap: anywhere; }
.card h3 { font-size: .95rem; margin: 0 0 6px; }
.card p { color: var(--muted); font-size: .8rem; margin: 2px 0; }
svg { width: 100%; height: auto; }
"""


def drawQuiver(data, size = 240):
    """One quiver as an inline SVG: vertices labelled, arrows with heads."""
    algebra = shapeKeys.deserialise(data)
    quiver = algebra.quiver
    graph = nx.Graph(quiver.to_undirected())
    if graph.number_of_nodes() > 1:
        raw = nx.spring_layout(graph, seed = 0)
    else:
        raw = {vertex: (0.0, 0.0) for vertex in graph.nodes}
    xs = [p[0] for p in raw.values()] or [0.0]
    ys = [p[1] for p in raw.values()] or [0.0]
    margin = 22
    span = max(max(xs) - min(xs), max(ys) - min(ys), 1e-9)
    scale = (size - 2 * margin) / span
    at = {v: (margin + (x - min(xs)) * scale, margin + (y - min(ys)) * scale)
          for v, (x, y) in raw.items()}
    parts = ['<svg viewBox="0 0 {0} {0}" xmlns="http://www.w3.org/2000/svg" role="img">'.format(size),
             '<defs><marker id="head" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" '
             'markerHeight="6" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="var(--line)"/>'
             '</marker></defs>']
    bundles = {}
    for tail, head, key in sorted(quiver.edges(keys = True)):
        bundles.setdefault((tail, head), []).append(key)
    radius = 9
    for (tail, head), keys in bundles.items():
        (x1, y1), (x2, y2) = at[tail], at[head]
        length = math.hypot(x2 - x1, y2 - y1) or 1.0
        ux, uy = (x2 - x1) / length, (y2 - y1) / length
        sx, sy = x1 + ux * radius, y1 + uy * radius
        ex, ey = x2 - ux * (radius + 2), y2 - uy * (radius + 2)
        for position, _key in enumerate(keys):
            bend = (position - (len(keys) - 1) / 2) * 18
            cx, cy = (sx + ex) / 2 - uy * bend, (sy + ey) / 2 + ux * bend
            parts.append('<path d="M{0:.1f},{1:.1f} Q{2:.1f},{3:.1f} {4:.1f},{5:.1f}" '
                         'fill="none" stroke="var(--line)" stroke-width="1.4" '
                         'marker-end="url(#head)"/>'.format(sx, sy, cx, cy, ex, ey))
    for vertex, (x, y) in at.items():
        parts.append('<circle cx="{0:.1f}" cy="{1:.1f}" r="{2}" fill="var(--node)" '
                     'stroke="var(--accent)" stroke-width="1.4"/>'.format(x, y, radius))
        parts.append('<text x="{0:.1f}" y="{1:.1f}" font-size="9" text-anchor="middle" '
                     'fill="var(--fg)">{2}</text>'.format(x, y + 3, html.escape(str(vertex))))
    parts.append('</svg>')
    return ''.join(parts)


def render(title, sections):
    """A whole page: a heading per section, a card per shape."""
    body = ['<h1>{0}</h1>'.format(html.escape(title))]
    for heading, entries in sections:
        body.append('<h2>{0}</h2><div class="grid">'.format(html.escape(heading)))
        for entry in entries:
            body.append('<div class="card"><h3>{0}</h3>{1}{2}</div>'.format(
                html.escape(entry['heading']), drawQuiver(entry['quiver']),
                ''.join('<p>{0}</p>'.format(html.escape(line)) for line in entry['lines'])))
        body.append('</div>')
    return ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width, initial-scale=1">'
            '<title>{0}</title><style>{1}</style></head><body>{2}</body></html>'
            .format(html.escape(title), _STYLE, ''.join(body)))


def sectionsFrom(tables, level, top):
    """Hubs, bridges and leftover hubs at one level, as page sections."""
    key = 'key{0}'.format(level)
    measures = shapeAtlas.shapeMeasures(tables, level).filter(~pl.col('isLine'))
    nodes = tables['nodes']

    def entries(frame):
        found = []
        for row in frame.head(top).iter_rows(named = True):
            quiver = json.loads(nodes.filter(pl.col(key) == row[key])['quiver'][0])
            found.append({
                'heading': '{0} classes, {1} starts'.format(row['classes'], row['starts']),
                'quiver': quiver,
                'lines': [shapeKeys.describe(shapeKeys.deserialise(quiver)),
                          'return rate {0:.2f}, first at depth {1}, leftover share {2:.2f}'.format(
                              row['returnRate'], row['firstDepth'], row['leftoverShare'])],
            })
        return found

    return [
        ('Hubs', entries(measures)),
        ('Bridges', entries(measures.filter(pl.col('bridge')))),
        ('Hubs among leftovers', entries(measures.filter(pl.col('leftoverShare') > 0)
                                         .sort('leftoverShare', 'starts', descending = True))),
    ]
