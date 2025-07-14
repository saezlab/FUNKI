# Building vars
VER = $(shell python3 -c "from funki import __version__; print(__version__.replace('.', '_'))")
ICON = assets/funki_favicon.ico
FLAGS = --onefile \
		--hidden-import=PIL._tkinter_finder \
		--hidden-import=scipy._cyutility \
		--recursive-copy-metadata scanpy \
		--clean \
		--distpath=./ \
		--name=FUNKI_v$(VER) \
		--icon=$(ICON)

build: compile clean

compile:
	pyinstaller ./src/funki/app/main.py $(FLAGS)

clean:
	rm -r *.spec build/
