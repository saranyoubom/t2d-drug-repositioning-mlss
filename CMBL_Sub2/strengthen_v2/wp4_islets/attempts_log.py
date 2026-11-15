"""attempts_log.py - append-only record of every download / access attempt of WP4 (protocol 4.3 and 4.6), written row by row so that an interruption loses nothing.
File: attempts.csv beside this script. Columns: when, dataset, item, url, size_bytes, label_source, outcome, reason. No side effects at import.
Usage: from attempts_log import record; record('GSE86469', 'series matrix', url, 38464, 'sample characteristics', 'downloaded', '')"""
import os, csv, time

PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'attempts.csv')
COLUMNS = ['when', 'dataset', 'item', 'url', 'size_bytes', 'label_source', 'outcome', 'reason']


def record(dataset, item, url, size_bytes, label_source, outcome, reason=''):
    new = not os.path.exists(PATH)
    with open(PATH, 'a', newline='', encoding='utf8') as f:
        w = csv.writer(f)
        if new: w.writerow(COLUMNS)
        w.writerow([time.strftime('%Y-%m-%d %H:%M:%S'), dataset, item, url, '' if size_bytes is None else size_bytes, label_source, outcome, reason])
