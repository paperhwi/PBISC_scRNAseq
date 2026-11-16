// Place the individually exported panels on a 178 x 258 mm figure page (template size of the PBISC figure set).
// Panel letters, row headings and notes are editable Arial text; panels are pictures (Change Picture with the same-name file).
const pptxgen = require('pptxgenjs'); const fs = require('fs'); const path = require('path');
// Optional: theme colours via apply_theme.js (PPTX_SKILL). Without it the deck keeps Office's stock theme; panels are unaffected.
const applyTheme = process.env.PPTX_SKILL ? require(process.env.PPTX_SKILL + '/scripts/apply_theme.js').applyTheme : async () => {};
const DIR = 'panels_out'; const items = JSON.parse(fs.readFileSync(path.join(DIR, 'manifest.json')));
const THEME = { name: 'PBISC figures', headFontFace: 'Arial', bodyFontFace: 'Arial',
  colors: { dk1: '000000', lt1: 'FFFFFF', dk2: '272727', lt2: 'F2F2F2', accent1: '0F4D92', accent2: 'B64342',
            accent3: '767676', accent4: '3775BA', accent5: '42949E', accent6: '9A4D8E', hlink: '0F4D92', folHlink: '9A4D8E' } };
const pres = new pptxgen();
pres.defineLayout({ name: 'FIG_178x258', width: 7.0087, height: 10.158 }); pres.layout = 'FIG_178x258';
pres.theme = { headFontFace: 'Arial', bodyFontFace: 'Arial' }; pres.title = 'PBISC Figure 7 ver5';
const C = pres.SchemeColor;
pres.defineSlideMaster({ title: 'FIGURE_PAGE', background: { color: 'FFFFFF' },
  objects: [{ placeholder: { options: { name: 'title', type: 'title', x: 0.15, y: 0.1, w: 6.71, h: 0.30,
    fontSize: 9, bold: true, color: C.text1, align: 'left', valign: 'top', margin: 0 }, text: '' } }] });
pres.addSection({ title: 'Main figures' });
const s = pres.addSlide({ masterName: 'FIGURE_PAGE', sectionTitle: 'Main figures' });
s.addText('Figure 7. Bulk-inferred single cells from new patients extend rare-population associations in ESKD COVID-19', { placeholder: 'title' });
for (const it of items) {
  s.addImage({ path: path.join(DIR, it.file), x: it.x, y: it.y, w: it.w, h: it.h,
    altText: it.file.replace(/\.png$/, '').replace(/_/g, ' '), objectName: 'panel_' + it.file.replace(/\.png$/, '') });
  if (it.label) s.addText(it.label, { x: it.x - 0.01, y: it.y - 0.03, w: 0.22, h: 0.2, fontSize: 10, bold: true, color: C.text1,
    margin: 0, valign: 'top', isTextBox: true, objectName: 'label_' + it.label });
}
const rows = [[3.91, 'B ASC dividing: the published severity association, detected with fewer Real patients (e\u2013g)'],
              [5.87, 'NK CD56bright: an association with death that emerges only after adding PBISC patients (h\u2013j)'],
              [7.83, 'HSPC: the published severity association, recovered with PBISC patients (k\u2013m)']];
for (const [y, t] of rows) {
  s.addText(t, { x: 0.32, y, w: 6.4, h: 0.17, fontSize: 7.5, bold: true, color: C.text1, margin: 0, valign: 'top', isTextBox: true, objectName: 'heading_' + y });
}
if (fs.existsSync('legend.txt')) s.addNotes(fs.readFileSync('legend.txt', 'utf8'));
(async () => { const out = process.argv[2]; await pres.writeFile({ fileName: out }); await applyTheme(out, THEME); console.log('wrote', out); })();
