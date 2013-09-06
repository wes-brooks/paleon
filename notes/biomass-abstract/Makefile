biomass-abstract.pdf: biomass-abstract.md
	pandoc -N --template=mytemplate.tex --variable mainfont=Georgia --variable sansfont=Arial  --variable fontsize=12pt --variable version=1.11 biomass-abstract.md --latex-engine=xelatex --toc -o biomass-abstract.pdf


out2.pdf:
	pandoc -s -w context test.md -o out2.tex
	context out2.tex

out.pdf:
	pandoc test.md -o out.pdf


