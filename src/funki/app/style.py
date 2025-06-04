from tkinter import ttk

from funki import _colors

def load_style():

    style = ttk.Style()

    style.configure(
        '.',
        font='Arial',
    )
    style.configure(
        'TFrame',
        background=_colors['white'],
    )
#    style.configure(
#        'TNotebook',
#        background=_colors['white'],
#    )
#    style.configure(
#        'TLabel',
#        background=_colors['white'],
#    )

    return style

