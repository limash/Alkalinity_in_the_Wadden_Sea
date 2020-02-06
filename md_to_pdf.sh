pandoc --pdf-engine=xelatex  -V geometry:"top=2cm, bottom=1.5cm, left=2cm, right=2cm"\
                             -V colorlinks -V urlcolor=NavyBlue -V toccolor=Red\
                             -V fontsize=12pt\
                              altogether.md -o test.pdf