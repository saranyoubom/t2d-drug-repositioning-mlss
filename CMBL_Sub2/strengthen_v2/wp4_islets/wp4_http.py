"""wp4_http.py - polite plain-HTTPS helpers for WP4 (stdlib only, no side effects at import).
Rules of the work package: plain public HTTPS from NCBI / EBI with a plain tool User-Agent, at most 3 requests per second, nothing sent about the user, 403/404 and challenge pages are
reported and NEVER worked around, a file above 300 MB is not fetched without an explicit decision (head() gives the size first), downloads are untrusted data written to the folder given.
Functions: head(url) -> dict(status, size, type); get_text(url) -> str; list_index(url) -> [(name, href)] from an HTML directory index; download(url, path, max_mb=300) -> bytes written.
Every call returns or raises; a result that looks like an HTML challenge page for a data file is flagged (is_challenge) so the caller can record it as blocked."""
import os, re, sys, time, urllib.request, urllib.error

USER_AGENT = 'geo-benchmark/1.0'
MIN_GAP = 0.34
_LAST = [0.0]


def throttle():
    wait = _LAST[0] + MIN_GAP - time.time()
    if wait > 0: time.sleep(wait)
    _LAST[0] = time.time()


def _open(url, method='GET', timeout=120):
    throttle()
    return urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': USER_AGENT}, method=method), timeout=timeout)


def head(url, tries=3):
    """status, size in bytes (None if absent), content type; errors are returned as status strings, not raised"""
    for k in range(tries):
        try:
            r = _open(url, 'HEAD'); cl = r.headers.get('Content-Length'); out = dict(status=r.status, size=int(cl) if cl else None, type=r.headers.get('Content-Type'), url=r.geturl()); r.close(); return out
        except urllib.error.HTTPError as e:
            if e.code in (403, 404, 405, 400): return dict(status=e.code, size=None, type=None, url=url)
            time.sleep(3 + 3 * k)
        except Exception as e:
            last = repr(e); time.sleep(3 + 3 * k)
    return dict(status='error: ' + last, size=None, type=None, url=url)


def get_text(url, tries=3, timeout=120, maxbytes=20_000_000):
    for k in range(tries):
        try:
            r = _open(url, 'GET', timeout); data = r.read(maxbytes + 1); r.close()
            if len(data) > maxbytes: raise RuntimeError('response larger than %d bytes: %s' % (maxbytes, url))
            return data.decode('utf8', 'replace')
        except urllib.error.HTTPError as e:
            if e.code in (403, 404): raise
            time.sleep(3 + 3 * k)
        except RuntimeError: raise
        except Exception: time.sleep(3 + 3 * k)
    raise RuntimeError('request failed ' + url)


def list_index(url):
    """names and hrefs of an HTML directory index (GEO ftp over https, EBI ftp)"""
    html = get_text(url)
    return [(re.sub(r'<[^>]+>', '', t).strip(), h) for h, t in re.findall(r'<a href="([^"?][^"]*)">(.*?)</a>', html, flags=re.S) if not h.startswith(('/', 'http', '?'))]


def is_challenge(first_bytes):
    s = first_bytes[:2000].lower()
    return b'<html' in s or b'<!doctype html' in s or b'recaptcha' in s


def download(url, path, max_mb=300, tries=3, chunk=1 << 20):
    """stream url into path (atomic .part rename). Refuses (RuntimeError) when the announced size is above max_mb; raises on 403/404; the first chunk of a data file that looks like an HTML page is
    saved as path + '.html_challenge' and reported as blocked (not worked around). Returns bytes written."""
    if os.path.exists(path) and os.path.getsize(path) > 0: return os.path.getsize(path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    h = head(url)
    if isinstance(h['status'], int) and h['status'] >= 400: raise urllib.error.HTTPError(url, h['status'], 'HEAD status %s' % h['status'], None, None)
    if h['size'] is not None and h['size'] > max_mb * 1e6: raise RuntimeError('size %.1f MB above the limit of %d MB: %s' % (h['size'] / 1e6, max_mb, url))
    last = None
    for k in range(tries):
        try:
            r = _open(url, 'GET', 180); n = 0; first = True
            with open(path + '.part', 'wb') as f:
                while True:
                    b = r.read(chunk)
                    if not b: break
                    if first:
                        first = False
                        if is_challenge(b) and not path.endswith(('.html', '.htm')):
                            open(path + '.html_challenge', 'wb').write(b[:5000]); r.close(); os.remove(path + '.part'); raise PermissionError('HTML challenge page instead of data: ' + url)
                    f.write(b); n += len(b)
                    if n > max_mb * 1e6 * 1.05: raise RuntimeError('download exceeded the size limit: ' + url)
            r.close(); os.replace(path + '.part', path); return n
        except urllib.error.HTTPError as e:
            last = e
            if e.code in (403, 404): raise
            time.sleep(3 + 3 * k)
        except (PermissionError, RuntimeError): raise
        except Exception as e:
            last = e; time.sleep(3 + 3 * k)
    raise RuntimeError('download failed %s (%r)' % (url, last))


if __name__ == '__main__':
    sys.stdout.reconfigure(encoding='utf8')
    if len(sys.argv) > 2 and sys.argv[1] == 'head': print(head(sys.argv[2]))
    elif len(sys.argv) > 2 and sys.argv[1] == 'index':
        for n, h in list_index(sys.argv[2]): print(n, '|', h)
