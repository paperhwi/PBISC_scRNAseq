import re
import v4_segments as V
src=open("build_deck_v3.py",encoding="utf-8").read()
parts=re.split(r"(?m)^(?=# ========================= )",src)
head=parts[0]; segs={p.split("\n")[0].replace("# ========================= ",""):p for p in parts[1:]}
order=list(segs.keys()); print(order)
tail_split=segs["24 discussion"].split("os.makedirs",1); save="os.makedirs"+tail_split[1]
out=[head]
for k in order:
    if k=="1 title": out.append(V.TITLE); continue
    if k=="2 summary": out.append(V.SUMMARY); continue
    if k=="21 synthesis": out.append(V.SYN); continue
    if k=="22 actions": out.append(V.ACT); continue
    if k=="23 limitations": out.append(V.LIM); continue
    if k=="24 discussion": out.append(V.DISC); out.append(V.APPX); continue
    out.append(segs[k])
    if k=="4 methods": out.append(V.PROV)
    if k=="9 volcano": out.append(V.LYMPH); out.append(V.MATCHED)
    if k=="15 caution IEG": out.append(V.DC)
    if k=="20b B1 specificity law": out+= [V.MYE_S,V.MYE_C,V.MK1,V.MK2]
save=save.replace("v3.pptx","v4.pptx").replace("speaker_notes_v3.md","speaker_notes_v4.md").replace("(v3)","(v4)")
save=save.replace("2026년 10월 27일 · {len(NOTE_LOG)}슬라이드 · 목표 40분 (질의 제외 약 37분)","2026년 10월 30일 · {len(NOTE_LOG)}슬라이드 · 본편 약 45분 + 부록 (질의응답용)")
out.append(save)
code="".join(out)
open("build_deck_v4.py","w",encoding="utf-8").write(code)
print("written", code.count("new()"))
