"""fda_fetch.py - cached, rate-limited, retrying client for the openFDA drug adverse event endpoint (FAERS).
Public data, plain HTTPS GET, no key, nothing about the user is sent (generic User-Agent). At most 2.5 requests per second;
openFDA allows 240 per minute and 1000 per day without a key, so a long run is split over days by simply re-running it
(every response is cached on disk and never requested twice).
Every response is stored as JSON in results/raw/ (file name = sha1 of the URL) together with its URL and retrieval time,
so that the analysis reads only saved files and can be re-run offline. HTTP 404 is what openFDA answers when a search
matches nothing; it is stored as a total of 0.
Usage as a module: from fda_fetch import total, counts, set_raw_dir ; python -I fda_fetch.py probe (syntax probe, writes to results/probe_raw/)."""
import os, sys, json, time, hashlib, urllib.request, urllib.parse, urllib.error

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = 'https://api.fda.gov/drug/event.json'
MIN_GAP = 0.4
_state = {'last': 0.0, 'raw': os.path.join(HERE, 'results', 'raw'), 'requests': 0, 'cached': 0, 'offline': False}


def set_raw_dir(path):
    _state['raw'] = path


def set_offline(flag=True):
    """offline: answer only from saved files and raise if a response is missing (used by the analysis)"""
    _state['offline'] = bool(flag)


def url_for(search=None, count=None, limit=1):
    parts = []
    if search:
        parts.append('search=' + urllib.parse.quote_plus(search))
    if count:
        parts.append('count=' + count)
    parts.append('limit=%d' % limit)
    return BASE + '?' + '&'.join(parts)


def _save(path, url, resp):
    rec = {'url': url, 'retrieved': time.strftime('%Y-%m-%d %H:%M:%S %z'), 'response': resp}
    tmp = path + '.tmp'
    with open(tmp, 'w', encoding='utf8') as f:
        json.dump(rec, f)
    os.replace(tmp, path)


def fetch(url, tries=9):
    raw = _state['raw']
    os.makedirs(raw, exist_ok=True)
    path = os.path.join(raw, hashlib.sha1(url.encode()).hexdigest() + '.json')
    if os.path.exists(path):
        _state['cached'] += 1
        with open(path, encoding='utf8') as f:
            return json.load(f)['response']
    if _state['offline']:
        raise RuntimeError('offline and not saved: ' + url)
    delay = 3.0
    last_err = None
    for k in range(tries):
        gap = time.time() - _state['last']
        if gap < MIN_GAP:
            time.sleep(MIN_GAP - gap)
        _state['last'] = time.time()
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'research-screen/1.0'})
            with urllib.request.urlopen(req, timeout=120) as r:
                resp = json.load(r)
            _state['requests'] += 1
            _save(path, url, resp)
            return resp
        except urllib.error.HTTPError as e:
            if e.code == 404:
                resp = {'meta': {'results': {'total': 0}}, 'results': [], 'note': 'HTTP 404: no match'}
                _state['requests'] += 1
                _save(path, url, resp)
                return resp
            last_err = 'HTTP %s' % e.code
            if e.code not in (429, 500, 502, 503, 504):
                raise
        except Exception as e:  # time-outs, DNS and SSL errors
            last_err = '%s: %s' % (type(e).__name__, str(e)[:120])
        print('  retry %d/%d after %s (wait %.0f s)' % (k + 1, tries, last_err, delay), flush=True)
        time.sleep(delay)
        delay = min(delay * 2, 120)
    raise RuntimeError('giving up on %s: %s' % (url, last_err))


def total(search=None):
    """number of reports matching the search (all reports when search is None)"""
    return fetch(url_for(search, None, 1))['meta']['results']['total']


def last_updated():
    return fetch(url_for(None, None, 1))['meta'].get('last_updated')


def counts(search, field, limit=1000):
    """top terms of a field among the reports matching the search: list of {term, count}"""
    return fetch(url_for(search, field, limit))['results']


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'probe':
        set_raw_dir(os.path.join(HERE, 'results', 'probe_raw'))
        D = 'patient.drug.openfda.substance_name'
        R = 'patient.reaction.reactionmeddrapt.exact'
        tests = [
            ('all reports', None),
            ('metformin', D + ':"metformin"'),
            ('NAUSEA exact', R + ':"NAUSEA"'),
            ('metformin AND NAUSEA exact', D + ':"metformin" AND ' + R + ':"NAUSEA"'),
            ('OR group exact', R + ':("NAUSEA" OR "VOMITING")'),
            ('OR group with a two-word PT', R + ':("NAUSEA" OR "DECREASED APPETITE")'),
            ('metformin AND OR group', D + ':"metformin" AND ' + R + ':("NAUSEA" OR "VOMITING")'),
            ('two drugs (metformin, sitagliptin)', D + ':"metformin" AND ' + D + ':"sitagliptin"'),
            ('two drugs AND OR group', D + ':"metformin" AND ' + D + ':"sitagliptin" AND ' + R + ':("NAUSEA" OR "VOMITING")'),
        ]
        for name, s in tests:
            t0 = time.time()
            try:
                print('%-40s %s   (%.1f s)' % (name, total(s), time.time() - t0), flush=True)
            except Exception as e:
                print('%-40s FAILED %s' % (name, e), flush=True)
        print('last_updated', last_updated())
