"""Tiny HTTP server with HTTP Range support, for streaming the rendered film.

Range support lets a browser seek the mp4 and start playback before the whole
file downloads. Serves animations/out/ on the chosen host:port.

    python3 serve.py --host 100.102.222.15 --port 8723
"""
import argparse
import http.server
import os
import re
import socketserver

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


class RangeHandler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *a, **k):
        super().__init__(*a, directory=ROOT, **k)

    def do_GET(self):
        path = self.translate_path(self.path)
        if os.path.isdir(path):
            return super().do_GET()
        try:
            f = open(path, "rb")
        except OSError:
            self.send_error(404)
            return
        try:
            size = os.fstat(f.fileno()).st_size
            ctype = self.guess_type(path)
            rng = self.headers.get("Range")
            if not rng:
                self.send_response(200)
                self.send_header("Content-Type", ctype)
                self.send_header("Content-Length", str(size))
                self.send_header("Accept-Ranges", "bytes")
                self.end_headers()
                self._copy(f, 0, size)
                return
            m = re.match(r"bytes=(\d+)-(\d*)", rng)
            start = int(m.group(1))
            end = int(m.group(2)) if m.group(2) else size - 1
            end = min(end, size - 1)
            length = end - start + 1
            self.send_response(206)
            self.send_header("Content-Type", ctype)
            self.send_header("Accept-Ranges", "bytes")
            self.send_header("Content-Range", f"bytes {start}-{end}/{size}")
            self.send_header("Content-Length", str(length))
            self.end_headers()
            self._copy(f, start, length)
        finally:
            f.close()

    def _copy(self, f, start, length):
        f.seek(start)
        remaining = length
        while remaining > 0:
            chunk = f.read(min(65536, remaining))
            if not chunk:
                break
            try:
                self.wfile.write(chunk)
            except (BrokenPipeError, ConnectionResetError):
                break
            remaining -= len(chunk)


class Server(socketserver.ThreadingTCPServer):
    allow_reuse_address = True
    daemon_threads = True


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", default="0.0.0.0")
    ap.add_argument("--port", type=int, default=8723)
    args = ap.parse_args()
    with Server((args.host, args.port), RangeHandler) as s:
        print(f"serving {ROOT} on http://{args.host}:{args.port}", flush=True)
        s.serve_forever()
