"""病原监测环状条码拼接 HTTP 服务（仅依赖标准库）。

环境变量 API_PORT 指定监听端口（默认 8080）。
路由：
  GET  /health   -> {"status": "ok"}
  POST /assemble -> {"sequences": [...]} -> unique / ambiguous / no_solution
"""
from __future__ import annotations

import json
import os
import traceback
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

from barcode import ValidationError, assemble

API_PORT = int(os.environ.get("API_PORT", "8080"))


class Handler(BaseHTTPRequestHandler):
    server_version = "BarcodeAPI/1.0"

    def _send_json(self, code: int, payload: dict) -> None:
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):  # noqa: N802
        path = self.path.split("?", 1)[0]
        if path == "/health":
            self._send_json(200, {"status": "ok", "service": "barcode-assembler"})
        else:
            self._send_json(404, {"error": "not_found", "path": path})

    def do_POST(self):  # noqa: N802
        path = self.path.split("?", 1)[0]
        if path != "/assemble":
            self._send_json(404, {"error": "not_found", "path": path})
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
            raw = self.rfile.read(length) if length else b""
            data = json.loads(raw.decode("utf-8")) if raw else {}
            if not isinstance(data, dict) or "sequences" not in data:
                raise ValidationError("请求体必须是包含 sequences 字段的 JSON 对象")
            result = assemble(data["sequences"])
            self._send_json(200, result)
        except ValidationError as exc:
            self._send_json(400, {"status": "invalid_input", "error": str(exc)})
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            self._send_json(400, {"status": "invalid_input", "error": f"请求体不是合法 JSON: {exc}"})
        except Exception:  # noqa: BLE001
            traceback.print_exc()
            self._send_json(500, {"status": "internal_error", "error": "服务内部错误"})

    def log_message(self, fmt, *args):  # 安静一点
        pass


def main() -> None:
    server = ThreadingHTTPServer(("0.0.0.0", API_PORT), Handler)
    print(f"barcode-assembler listening on 0.0.0.0:{API_PORT}", flush=True)
    server.serve_forever()


if __name__ == "__main__":
    main()
