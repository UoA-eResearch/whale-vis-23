import numpy as np
from bokeh.events import DocumentReady
from bokeh.io import curdoc
from bokeh.models import PolyAnnotation, Label, Plot, BoxAnnotation, CustomJS
from datetime import datetime


def date_annotation(figure: Plot, date: datetime, x_pos=0.9, y_pos=0.95, text_size=20):
    """Add a date annotation to a map figure"""
    # Place on the figure
    abs_x_pos = x_pos * figure.width
    abs_y_pos = y_pos * figure.height

    # Label the date
    date_label = Label(x=abs_x_pos, y=abs_y_pos,
                        text=date.strftime('%d-%m-%Y'), text_font_size=f'{text_size}pt',
                        text_align='right', text_baseline='top',
                        x_units='screen', y_units='screen')
    figure.add_layout(date_label)


def north_arrow(figure: Plot, x_pos=0.05, y_pos=0.95, arrow_size=30, text_size=20):
    """Add a north arrow to a map figure"""
    # North arrow
    xs = np.array([0, .5, 1, .5])
    ys = np.array([0, 1, 0, .5])
    # Center
    xs -= 0.5
    # Scale
    xs *= arrow_size
    ys *= arrow_size
    # Place on the figure
    abs_x_pos = x_pos * figure.width
    abs_y_pos = y_pos * figure.height - arrow_size
    xs += abs_x_pos
    ys += abs_y_pos

    northArrow = PolyAnnotation(fill_color='black', fill_alpha=1,
                                xs=list(xs), ys=list(ys),
                                xs_units='screen', ys_units='screen')
    northLetter = Label(x=abs_x_pos, y=abs_y_pos,
                        text='N', text_font_size=f'{text_size}pt',
                        text_align='center', text_baseline='top',
                        x_units='screen', y_units='screen')
    figure.add_layout(northArrow)
    figure.add_layout(northLetter)


def scale_bar(figure: Plot, halved=False, convert_from_deg=False):
    """Add a scale bar to a map figure"""
    # Style either half black/white, or solid black
    # Recommend not using halved because the line needed around the outside edge is on the outside of the bar, making it
    # bigger than it should be
    if halved:
        scalebar_h0 = BoxAnnotation(left=0, right=1, top=0, bottom=1, fill_color='black', fill_alpha=0, line_color='black',
                                    line_width=2, line_alpha=1)
        scalebar_h1 = BoxAnnotation(left=0, right=1, top=0, bottom=1, fill_color='black', fill_alpha=1, line_color='black',
                                    line_width=2, line_alpha=1)
    else:
        scalebar_h0 = BoxAnnotation(left=0, right=1, top=0, bottom=1, fill_color='black', fill_alpha=1, line_color='black',
                                    line_width=2, line_alpha=0)
        scalebar_h1 = BoxAnnotation(left=0, right=1, top=0, bottom=1, fill_color='black', fill_alpha=1, line_color='black',
                                    line_width=2, line_alpha=0)

    # Label the bar length
    scalebar_l = Label(x=0, y=0, text='100km', text_font_size='14pt', text_align='right', text_baseline='top')
    figure.add_layout(scalebar_h0)
    figure.add_layout(scalebar_h1)
    figure.add_layout(scalebar_l)

    # JS code that calculates an appropriate length for the bar (2 s.f.), places it
    # then places the label and updates it's text
    slider_callback = CustomJS(
        args={'sb': [scalebar_h0, scalebar_h1], 'sl': scalebar_l, 'x': figure.x_range, 'y': figure.y_range}, code=f'''
    const lerp = (a, b, v) => a + v * (b - a);
    const [x0, x1, y0, y1] = [x.start, x.end, y.start, y.end]

    let width = (0.2 * (x1-x0)).toPrecision(2)

    function haversine_distance(lat1, lon1, lat2, lon2) {{
        const R = 6371e3; // metres
        const la1 = lat1 * Math.PI/180; // convert to radians
        const la2 = lat2 * Math.PI/180;
        const dlat = (lat2-lat1) * Math.PI/180;
        const dlon = (lon2-lon1) * Math.PI/180;

        const a = Math.sin(dlat/2) * Math.sin(dlat/2) +
                  Math.cos(la1) * Math.cos(la2) *
                  Math.sin(dlon/2) * Math.sin(dlon/2);
        const c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1-a));

        return R * c;
    }}
    
    let width_text = {'haversine_distance(y0, x1 - width, y0, x1).toPrecision(2)' if convert_from_deg else 'width'};

    let right = lerp(x0, x1, 0.95);
    let left = right - width;
    let middle = lerp(left, right, 0.5);
    let top = lerp(y0, y1, 0.05);
    let bottom = lerp(y0, y1, 0.04);

    sb[0].left = left;
    sb[0].right = middle;
    sb[0].top = top;
    sb[0].bottom = bottom;

    sb[1].left = middle;
    sb[1].right = right;
    sb[1].top = top;
    sb[1].bottom = bottom;

    sl.x = right;
    sl.y = bottom;
    if (width_text > 10000) {{
        sl.text = `${{width_text / 1000}} km`
    }} else {{
        sl.text = `${{Number(width_text)}} m`
    }}
    ''')

    figure.x_range.js_on_change('start', slider_callback)
    figure.x_range.js_on_change('end', slider_callback)
    figure.y_range.js_on_change('start', slider_callback)
    figure.y_range.js_on_change('end', slider_callback)

    curdoc().js_on_event(DocumentReady, slider_callback)
