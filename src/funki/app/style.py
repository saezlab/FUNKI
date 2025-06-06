from tkinter import ttk

from funki import _colors

def load_style(root):

    style = ttk.Style(root)

    style.configure(
        '.',
        font=('Arial', 12),
    )
    # Tabs
    style.configure(
        'TNotebook',
        tabposition='wn',
    )
    style.configure(
        'TNotebook.Tab',
        padding=(5, 5, 5, 5),
        background=_colors['blue'],
        foreground=_colors['white'],
        width=12,
        anchor='EW',
    )
    style.map(
        'TNotebook.Tab',
        background=[
            ('selected', _colors['teal']),
            ('active', _colors['aqua']),
            #('disabled', _colors['teal'])
        ],
        foreground=[
            ('selected', _colors['yellow']),
            ('active', _colors['white']),
            #('disabled', _colors['yellow'])
        ],
    )

    return style


table = ttk.Style()

table.configure(
    'cell',
    anchor='right',
    background=_colors['white'],
    font=('Courier', 10),
)
table.configure(
    'index',
    anchor='left',
    background=_colors['gray'],
    font=('Courier', 10, 'bold'),
)
table.configure(
    'column',
    anchor='center',
    background=_colors['gray'],
    font=('Courier', 10, 'bold'),
)