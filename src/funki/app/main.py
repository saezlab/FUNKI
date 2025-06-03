import tkinter as tk

from funki import __version__

class Funki(tk.Frame):
    def __init__(self, root=None, **options):
        super().__init__(root)
        self.root = root
        self.root.title('FUNKI v%s' % __version__)


if __name__ == '__main__':
    root = tk.Tk()

    app = Funki(root=root)
    app.mainloop()