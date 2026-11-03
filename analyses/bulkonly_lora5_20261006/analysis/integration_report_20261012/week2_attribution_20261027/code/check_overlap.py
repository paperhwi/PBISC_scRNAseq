"""Geometry audit: shapes outside the slide and text/picture/table collisions."""
from pptx import Presentation
from pptx.util import Emu
P=Presentation("deck/PBISC 질병 DEG 희귀집단 검증 보고서 v2.pptx")
W,H=P.slide_width,P.slide_height
E=lambda v: v/914400
def rect(sh): return (sh.left,sh.top,sh.left+sh.width,sh.top+sh.height)
def inter(a,b):
    x=min(a[2],b[2])-max(a[0],b[0]); y=min(a[3],b[3])-max(a[1],b[1])
    return (x,y) if x>0 and y>0 else None
bad=0
for i,s in enumerate(P.slides,1):
    items=[]
    for sh in s.shapes:
        if sh.left is None: continue
        r=rect(sh)
        if r[0]<-9144 or r[1]<-9144 or r[2]>W+9144 or r[3]>H+9144:
            print(f"S{i}: OUT OF BOUNDS {sh.shape_type} {[round(E(v),2) for v in r]}"); bad+=1
        kind=("pic" if sh.shape_type==13 else "tbl" if sh.has_table else
              "txt" if (sh.has_text_frame and sh.text_frame.text.strip()) else "box")
        items.append((kind,r,(sh.text_frame.text[:22] if sh.has_text_frame else "")))
    cont=[x for x in items if x[0] in ("pic","tbl","txt")]
    for a in range(len(cont)):
        for b in range(a+1,len(cont)):
            ka,ra,ta=cont[a]; kb,rb,tb=cont[b]
            ov=inter(ra,rb)
            if ov and ov[0]>Emu(0.05*914400) and ov[1]>Emu(0.05*914400):
                print(f"S{i}: {ka}({ta!r}) x {kb}({tb!r}) overlap {E(ov[0]):.2f}x{E(ov[1]):.2f} in")
                bad+=1
print("ISSUES:",bad)
