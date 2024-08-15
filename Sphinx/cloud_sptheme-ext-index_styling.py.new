"""cloud_sptheme.ext.index_styling - improves css styling for genindex"""
import logging; log = logging.getLogger(__name__)
import re
#jkfrom jinja2 import Markup as literal, escape
from jinja2.utils import markupsafe as literal
from cloud_sptheme import __version__

prefix = r"^(?P<name>.*)\("
suffix = r"\)$"
_attr_re = re.compile(prefix + r"(?P<left>)(?P<loc>.*)(?P<right> attribute)" + suffix)
_meth_re = re.compile(prefix + r"(?P<left>)(?P<loc>.*)(?P<right> method)" + suffix)
_fc_re = re.compile(prefix + r"(?P<left>class in |in module )(?P<loc>.*)(?P<right>)" + suffix)
_mod_re = re.compile(prefix + r"module" + suffix)

def format_index_name(name):
    while True:
        m = _attr_re.match(name)
        if m:
            name, left, loc, right = m.group("name","left", "loc", "right")
            type = "attribute"
            break
        m = _meth_re.match(name)
        if m:
            name, left, loc, right = m.group("name","left", "loc", "right")
            type = "method"
            break
        m = _fc_re.match(name)
        if m:
            name, left, loc, right = m.group("name","left", "loc", "right")
            if left.startswith("class"):
                type = "class"
            else:
                type = "function"
            break
        m = _mod_re.match(name)
        if m:
            name = m.group("name")
            left = "module"
            loc = right = ''
            type = "module"
            break
        return name
    if loc:
        #jk loc = literal('<span class="location">') + escape(loc) + literal("</span>")
        loc = literal.Markup('<span class="location">') + literal.escape(loc) + literal.Markup("</span>")
    cat = left + loc + right
    return literal.escape(name) + literal.Markup('<span class="category ' + type + '">') + literal.escape(cat) + literal.Markup("</span>")

def mangle_index(app, pagename, templatename, ctx, event_arg):
    if pagename != "genindex":
        return
    fmt = format_index_name
    for key, entries in ctx['genindexentries']:
        for idx, entry in enumerate(entries):
            name, data = entry
            entries[idx] = fmt(name), data
            # NOTE: data is list of [links, subitems, some_key],
            #       though 'key' not added until sphinx 1.4 (1.3?)
            subitems = data[1]
            for idx, entry in enumerate(subitems):
                name, links = entry
                subitems[idx] = fmt(name), links

def setup(app):
    app.connect('html-page-context', mangle_index)

    # identifies the version of our extension
    return {'version': __version__}
