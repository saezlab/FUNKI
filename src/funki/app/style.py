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
    style.configure(
        'TNotebook.Tab',
        padding=(5, 5, 5, 5),
        background=_colors['red']
    )
    style.map(
        'TNotebook.Tab',
        background=[
            ('selected', _colors['blue']),
            ('active', _colors['teal']),
            ('disabled', _colors['teal'])
        ],
        foreground=[
            ('selected', _colors['red']),
            ('active', _colors['yellow']),
            ('disabled', _colors['yellow'])
        ],
    )
#    style.configure(
#        'TLabel',
#        background=_colors['white'],
#    )

    return style

