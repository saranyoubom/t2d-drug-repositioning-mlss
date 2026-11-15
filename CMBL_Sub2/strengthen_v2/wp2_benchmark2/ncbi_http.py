"""ncbi_http.py - tiny politeness layer for every NCBI request of WP2 (stdlib only, no side effects at import).
throttle() keeps consecutive requests at least MIN_GAP seconds apart (about 3 per second, the limit set for this work package); http_get() is a GET with the plain tool User-Agent used in
part 1 ('geo-benchmark/1.0'), retries on transient errors, and raises on 404 and 403 (never retried, never worked around). No e-mail address or key is sent."""
import time, urllib.request, urllib.error

MIN_GAP = 0.34
_LAST = [0.0]
USER_AGENT = 'geo-benchmark/1.0'


def throttle():
    wait = _LAST[0] + MIN_GAP - time.time()
    if wait > 0: time.sleep(wait)
    _LAST[0] = time.time()


def http_get(url, binary=True, tries=4, timeout=180, stream=False, method='GET'):
    """GET (or HEAD) a URL. stream=True returns the open response (caller closes it); otherwise the body as bytes (binary) or text. HTTPError 404/403 is raised at once."""
    for k in range(tries):
        throttle()
        try:
            r = urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': USER_AGENT}, method=method), timeout=timeout)
            if stream or method == 'HEAD': return r
            data = r.read(); r.close(); return data if binary else data.decode('utf8', 'replace')
        except urllib.error.HTTPError as e:
            if e.code in (404, 403): raise
            time.sleep(3 + 3 * k)
        except Exception: time.sleep(3 + 3 * k)
    raise RuntimeError('download failed ' + url)
