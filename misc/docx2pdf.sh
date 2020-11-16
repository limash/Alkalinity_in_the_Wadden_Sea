# pandoc --pdf-engine=xelatex --include-in-header properties.tex\
#                             --variable classoption=twocolumn\
#                             --variable papersize=a4paper\
#                             -s main.docx -o main.pdf

pandoc --pdf-engine=xelatex -V geometry:"top=2cm, bottom=1.5cm, left=2cm, right=2cm"\
                            -V colorlinks -V urlcolor=NavyBlue -V toccolor=Red\
                            -V fontsize=12pt\
                            -s main.docx -o main.pdf