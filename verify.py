import urllib.request

resp = urllib.request.urlopen("http://127.0.0.1:8050")
html = resp.read().decode()

has_hs = "chart-hs-img" in html
has_th = "chart-th-img" in html
has_b64 = "data:image/png;base64" in html

print(f"Has chart-hs-img: {has_hs}")
print(f"Has chart-th-img: {has_th}")
print(f"Has base64 image: {has_b64}")
print(f"Result: {'PASS' if has_hs and has_th and has_b64 else 'FAIL'}")
