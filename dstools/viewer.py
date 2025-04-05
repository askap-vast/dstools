import warnings
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from astropy.wcs import FITSFixedWarning
from matplotlib.path import Path as Polypath
from matplotlib.widgets import Button
from mpl_point_clicker import clicker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from dstools.imaging import Image

warnings.filterwarnings("ignore", category=FITSFixedWarning, append=True)


@dataclass
class Viewer:
    images: list[Image]

    def __post_init__(self):
        for image in self.images:
            setattr(self, f"{image.name}_image", image)

        self.active_image = self.images[0]
        self.axes = []

        self._setup()

        # Set up mask array
        image_shape = self.images[0].data.shape
        self.mask = np.full(image_shape, True)
        span = np.arange(image_shape[0])
        self.coords = np.dstack(np.meshgrid(span, span)).reshape(-1, 2)

        # Plot first image
        self.image_object = self.ax.imshow(
            self.active_image.data,
            norm=self.active_image.norm,
            cmap="plasma",
        )

        mask = np.ma.masked_where(self.mask, ~self.mask) * 1
        self.mask_plot = self.ax.imshow(
            mask,
            cmap="bwr_r",
            alpha=0.5,
        )

        self.add_norm_buttons()
        self.add_click_handler()

        plt.show()

    def add_norm_buttons(self):
        """Add buttons for switching color scales."""

        scales = [0.9, 0.95, 0.99, 0.995, 0.999, 0.9999]
        text_ax = inset_axes(
            self.ax,
            width="40%",
            height="3%",
            loc="lower left",
            bbox_to_anchor=(0.74, -0.05, 1, 1),
            bbox_transform=self.ax.transAxes,
        )
        text_ax.text(
            0,
            0,
            "Colormap Normalisation Scale",
            horizontalalignment="center",
            size=11,
        )

        text_ax.set_xticks([])
        text_ax.set_yticks([])
        text_ax.set_frame_on(False)

        for i, scale in enumerate(scales):
            xpos = i * 0.08
            button_ax = inset_axes(
                self.ax,
                width="7%",
                height="7%",
                loc="lower left",
                bbox_to_anchor=(0.5 + xpos, -0.15, 1, 1),
                bbox_transform=self.ax.transAxes,
            )

            toggle = Button(button_ax, f"{scale * 100:.10g}%")
            toggle.on_clicked(self._update_colorscale(scale * 100))
            setattr(self, f"scale{scale}_toggle", toggle)

    def _update_colorscale(self, percentile):
        data = self.image_object.get_array().data
        vmin = np.nanpercentile(data, 100 - percentile)
        vmax = np.nanpercentile(data, percentile)

        def _update(val):
            self.image_object.norm.vmin = vmin
            self.image_object.norm.vmax = vmax

            # Redraw the figure to ensure it updates
            self.fig.canvas.draw_idle()

            return

        return _update

    def _setup(self):
        self.fig = plt.figure(figsize=(8, 9))
        self.ax = self.fig.add_subplot(111, projection=self.images[0].wcs)
        self.axes.append(self.ax)

        # Remove coordinates / ticklabels
        for axis in self.ax.coords:
            axis.set_axislabel("")
            axis.set_ticklabel_visible(False)
            axis.set_ticks_visible(False)

        # Add button toggle for each image
        for i, image in enumerate(self.images):
            xpos = i * 0.15
            _, toggle = self.add_button(image, xpos=xpos)
            setattr(self, f"{image.name}_toggle", toggle)

        # Set figure borders
        self.fig.subplots_adjust(
            bottom=0.15,
            left=0.01,
            right=0.99,
            top=0.97,
        )

    def add_click_handler(self):
        self.clicker = clicker(
            self.ax,
            ["mask"],
            markers=["o"],
            colors=["r"],
            linestyle="-",
            markersize=3,
            disable_legend=True,
        )

        self.fig.canvas.mpl_connect(
            "key_press_event",
            lambda event: self._on_press(event),
        )

    def _update_images(self):
        mask = np.ma.masked_where(self.mask, ~self.mask)
        image = self.active_image.data

        self.image_object.set_data(image)
        self.image_object.norm = self.active_image.norm

        self.mask_plot.set_data(mask * 1)

        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()

        return

    def _on_press(self, event):
        # TODO: refactor handling of image refreshing and try to speed up?

        if event.key not in ["x", "c"]:
            return

        points = self.clicker.get_positions()["mask"]

        if len(points) > 2:
            # Generate polygon from selected points to mask model
            polygon = Polypath(points)
            inside_polygon = polygon.contains_points(self.coords).reshape(
                self.mask.shape
            )

            if event.key == "x":
                self.mask = np.logical_and(self.mask, ~inside_polygon)
            elif event.key == "c":
                self.mask[inside_polygon] = True

        self._update_images()
        self.clicker.clear_positions()

        return

    # def add_slider(self):
    #     slider_ax = inset_axes(
    #         self.ax,
    #         width="30%",
    #         height="10%",
    #         loc="lower left",
    #         bbox_to_anchor=(0.6, -0.15, 1, 1),
    #         bbox_transform=self.ax.transAxes,
    #     )

    #     self.slider = RangeSlider(
    #         slider_ax,
    #         "Threshold",
    #         np.nanmin(self.active_image.data),
    #         np.nanmax(self.active_image.data),
    #     )

    #     self.slider.on_changed(self._update_colorscale())

    #     return

    def add_button(self, image, xpos):
        button_ax = inset_axes(
            self.ax,
            width="10%",
            height="10%",
            loc="lower left",
            bbox_to_anchor=(xpos, -0.15, 1, 1),
            bbox_transform=self.ax.transAxes,
        )

        toggle = Button(button_ax, image.name)
        toggle.on_clicked(self._switch_image(image))

        return button_ax, toggle

    def _switch_image(self, image):
        def callback(event):
            self.active_image = image
            self._update_images()

            return

        return callback
