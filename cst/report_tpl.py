"""Шаблон отчёта по ДЗ (кафедра РЛ1) из присланного примера: титул с исполнителями К.Д. Соломонов, О.Е. Смирнов,
проверяющий, СОДЕРЖАНИЕ (поле TOC). Тело примера удаляется, тема заменяется.
  from report_tpl import open_template
  d = open_template("ТЕМА ПРОПИСНЫМИ")
"""
import os
from docx import Document
from docx.oxml.ns import qn

HERE = os.path.dirname(os.path.abspath(__file__))
TPL = os.path.join(HERE, "tpl_dz.docx")      # копия Downloads\Отчёт_ДЗ_ФАР_ФВ_Фарадея_Ka.dot, сохранённая Word как .docx


def open_template(tema, kurs=None):
    d = Document(TPL)
    body = d.element.body
    els = list(body)
    # тема: абзац, начинающийся с ФЕРРИТОВЫЙ
    for p in d.paragraphs:
        if p.text.strip().startswith("ФЕРРИТОВЫЙ"):
            p.runs[0].text = tema
            for r in p.runs[1:]:
                r.text = ""
            break
    if kurs:                                   # строка «ПО КУРСУ …»: прогоны «ПО КУРСУ», « », название
        for p in d.paragraphs:
            if p.text.strip().startswith("ПО КУРСУ"):
                p.runs[-1].text = kurs
                for r in p.runs[2:-1]:
                    r.text = ""
                break
    # тело примера: всё после конца поля оглавления (fldChar end) до sectPr
    i_toc = next(i for i, e in enumerate(els) if "СОДЕРЖАНИЕ" in "".join(t.text or "" for t in e.iter(qn("w:t"))))
    depth, started, i_end = 0, False, None              # конец внешнего поля TOC (внутри есть поля PAGEREF)
    for i in range(i_toc, len(els)):
        for fc in els[i].iter(qn("w:fldChar")):
            t = fc.get(qn("w:fldCharType"))
            if t == "begin":
                depth += 1; started = True
            elif t == "end":
                depth -= 1
        if started and depth == 0:
            i_end = i; break
    for e in els[i_end + 1:]:
        if e.tag != qn("w:sectPr"):
            body.remove(e)
    return d
