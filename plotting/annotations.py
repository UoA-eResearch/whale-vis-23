import numpy as np
from bokeh.events import DocumentReady
from bokeh.io import curdoc
from bokeh.models import PolyAnnotation, Label, Plot, BoxAnnotation, CustomJS, Image
from datetime import datetime
from PIL import Image


def _lerp(a, b, v):
    """Linear interpolation"""
    return a + v * (b - a)


def date_annotation(figure: Plot, date: datetime, bounds, x_pos=0.9, y_pos=0.9, text_size=20):
    """Add a date annotation to a map figure"""
    # Place on the figure
    abs_x_pos = _lerp(bounds[0], bounds[2], x_pos)
    abs_y_pos = _lerp(bounds[1], bounds[3], y_pos)

    # Label the date
    date_label = Label(x=abs_x_pos, y=abs_y_pos,
                        text=date.strftime('%d-%m-%Y'), text_font_size=f'{text_size}pt',
                        text_align='center', text_baseline='bottom',
                        x_units='data', y_units='data')
    figure.add_layout(date_label)


def _pack_image(image):
    """Pack an rgba image into a 2d array"""
    # See https://docs.bokeh.org/en/latest/docs/examples/topics/images/image_rgba.html
    out = np.empty((image.shape[0], image.shape[1]), dtype=np.uint32)
    view = out.view(dtype=np.uint8).reshape((image.shape[0], image.shape[1], 4))
    view[:, :, :] = image
    return out


def add_logo(figure: Plot, filename, bounds, x_pos=0.15, y_pos=0.95, height=0.1, anchor='top_left'):
    # Load image from disk
    im = Image.open(filename)
    image = np.array(im.getdata()).reshape((im.size[1], im.size[0], 4))
    packed_image = _pack_image(image)
    packed_image = np.flipud(packed_image)
    im.close()

    # Convert relative screen position to absolute data position
    abs_x_pos = _lerp(bounds[0], bounds[2], x_pos)
    abs_y_pos = _lerp(bounds[1], bounds[3], y_pos)

    # Scale image to desired height
    image_aspect = packed_image.shape[1] / packed_image.shape[0]
    data_height = height * (bounds[3] - bounds[1])
    data_width = image_aspect * data_height

    figure.image_rgba(image=[packed_image], x=abs_x_pos, y=abs_y_pos,
                      dw=data_width, dh=data_height,
                      dw_units='data', dh_units='data',
                      anchor=anchor)


def north_arrow(figure: Plot, bounds, x_pos=0.05, y_pos=0.95, arrow_size=.05, text_size=20):
    """Add a north arrow to a map figure"""
    # North arrow
    xs = np.array([0, .5, 1, .5])
    ys = np.array([0, 1, 0, .5])
    # Center
    xs -= 0.5
    # Scale
    scale = arrow_size * (bounds[2] - bounds[0])
    xs *= scale
    ys *= scale
    # Place on the figure
    abs_x_pos = _lerp(bounds[0], bounds[2], x_pos)
    abs_y_pos = _lerp(bounds[1], bounds[3], y_pos) - scale
    xs += abs_x_pos
    ys += abs_y_pos

    northArrow = PolyAnnotation(fill_color='black', fill_alpha=1,
                                xs=list(xs), ys=list(ys),
                                xs_units='data', ys_units='data')
    northLetter = Label(x=abs_x_pos, y=abs_y_pos,
                        text='N', text_font_size=f'{text_size}pt',
                        text_align='center', text_baseline='top',
                        x_units='data', y_units='data')
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

    let right = lerp(x0, x1, 0.98);
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
