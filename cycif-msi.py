import os
import yaml

from tifffile import TiffFile
from tifffile import imread
import dask.array as da
import zarr

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns

from skimage.util import img_as_uint
from skimage.util import img_as_float
from skimage.color import gray2rgb

from natsort import natsorted

from scipy.stats import f_oneway
from bioinfokit.analys import stat

annotate_points = False
density = True
###############################################################################

# path to output directory
save_dir = 'results'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

# read unfiltered single-cell data
cycif = pd.read_csv('cycif_downsample.csv')

# slice unfiltered cells from current sample
clean_cells = cycif[~cycif['cluster'].isnull()]

# slice cluster 0 cells (CD73 dim cells) from current sample
cd73_dim = clean_cells[clean_cells['cluster'] == 0]

# slice cluster 2 cells (CD73 bright cells) from current sample
cd73_bright = clean_cells[clean_cells['cluster'] == 2]

# read upsampled adenosine image
adeno_upsample = pd.read_csv('msi_upsample.csv')
adeno_upsample = adeno_upsample.to_numpy()
adeno_upsample = img_as_uint(adeno_upsample)

# read raster height, width

# ensure upsampled adenosine image is floating point, grayscale
img = img_as_float(adeno_upsample)

# save 2D-grayscale adenosine image as a unique variable
adeno_gray = img

# adjust adenosine image contrast settings
bottom_cutoff = 0.05  # np.min(img)
top_cutoff = 0.8
img = np.clip(img, bottom_cutoff, top_cutoff)

# rescale 0-1
img = (img - np.min(img)) / (np.max(img) - np.min(img))

# apply colormap to adenosine image
img_rgb = plt.cm.Greys(img)  # auto converts to RGBA
img_rgb = img_rgb[:, :, 0:3]  # drop alpha values
img_rgb *= 1.0  # modify intensity

# plot
fig = plt.figure(figsize=(7.0, 3.0))

# grid specifications
gs = plt.GridSpec(
    nrows=3, ncols=6, figure=fig,
    height_ratios=[1.0, 1.0, 0.5],
    width_ratios=[0.5, 1.0, 1.0, 1.0, 0.75, 0.75]
    )

ax_img_legend = fig.add_subplot(gs[0, 0])
ax_adeno_img = fig.add_subplot(gs[0:3, 1:4])
ax_adeno_plot = fig.add_subplot(gs[0:2, 4:6])

plt.subplots_adjust(
    left=0.03, bottom=0.13, right=0.97,
    top=0.87, wspace=1.0, hspace=0.2
    )

# add merged images to figure
for ax, img, label in zip(
  [ax_adeno_img],
  [img_rgb],
  ['adenosine']
  ):

    im = ax.imshow(img, cmap='Greys', aspect='auto')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('bottom', size='3%', pad=0.0)
    cax.tick_params(labelsize=5)
    cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
    cbar.set_label('adenosine', fontsize=7)
    ax.axis('off')

# add image legend
legend_elements = []

legend_elements.append(
    Line2D([0], [0], marker='o',
           color='none',
           label='$\mathregular{CD73^{dim}}$ cells (cluster 0)',
           markerfacecolor='tab:blue',
           markeredgecolor='k',
           lw=0.005, markersize=7))

legend_elements.append(
    Line2D([0], [0], marker='o',
           color='none',
           label='$\mathregular{CD73^{bright}}$ cells (cluster 2)',
           markerfacecolor='tab:red',
           markeredgecolor='k',
           lw=0.005, markersize=7))

if annotate_points:
    legend_elements.append(
        Line2D([0], [0], marker='x',
               color='none',
               label='Quality grid box',
               markerfacecolor='lime',
               markeredgecolor='lime',
               lw=0.005, markersize=5))

ax_img_legend.legend(
    handles=legend_elements, prop={'size': 4},
    loc='upper left', bbox_to_anchor=[0.0, 0.98],
    labelspacing=1.5, frameon=False
    )

ax_img_legend.axis('off')

# create box grid
# get evenly spaced numbers over x and y axes ranges
# according to pixel-size of the original MSI images
raster_height = 30
raster_width = 36
num_height_samples = round(img_rgb.shape[0]/raster_height)
num_width_samples = round(img_rgb.shape[1]/raster_width)
row_rnge = np.linspace(
    0, img_rgb.shape[0], num_height_samples).astype(int)
col_rnge = np.linspace(
    0, img_rgb.shape[1], num_width_samples).astype(int)

# initialize dataframe to store per box measurements
overall_measurements = pd.DataFrame(
    {'row': pd.Series(dtype='int'),
     'column': pd.Series(dtype='int'),
     'cd73_bright': pd.Series(dtype='float'),
     'cd73_dim': pd.Series(dtype='float'),
     'concentration_adeno': pd.Series(dtype='float')}
     )

# partition image data into boxes and evaluate
row_start = 0
for e, row_end in enumerate(row_rnge[1:]):
    for ax in [ax_adeno_img]:
        ax.plot(
            [0, img_rgb.shape[1]], [row_end, row_end],
            c='grey', lw=0.04
            )
    col_start = 0
    for j, col_end in enumerate(col_rnge[1:]):
        for ax in [ax_adeno_img]:
            ax.plot(
                [col_end, col_end], [0, img_rgb.shape[0]],
                c='grey', lw=0.04
                )

        # get CD73 bright cells from current box
        box_cd73_bright = (
            cd73_bright[
                (cd73_bright['Y_centroid'].between(
                 row_start, row_end)) &
                (cd73_bright['X_centroid'].between(
                 col_start, col_end))]
            .copy()
            )

        # get CD73 dim cells from current box
        box_cd73_dim = (
            cd73_dim[
                (cd73_dim['Y_centroid'].between(
                 row_start, row_end)) &
                (cd73_dim['X_centroid'].between(
                 col_start, col_end))]
            .copy()
            )

        # get filtered cells from current box
        box_clean = (
            clean_cells[
                (clean_cells['Y_centroid'].between(
                 row_start, row_end))
                &
                (clean_cells['X_centroid'].between(
                 col_start, col_end))].copy()
            )

        # get total cells from current box
        box_total = (
            cycif[
                (cycif['Y_centroid'].between(
                 row_start, row_end))
                &
                (cycif['X_centroid'].between(
                 col_start, col_end))].copy()
            )

        # slice adenosine image to box size
        box_adeno = adeno_gray[row_start:row_end, col_start:col_end]

        # filter boxes for quality:
        # 1) apply lower cutoff on cell number per box
        # 2) apply fraction of filtered data per box
        # 3) remove boxes with NANs in the ratiometric images
        # caused by 0/0 errors. This also helps remove
        # boxes at the tissue boundaries.
        if (
            len(box_total) >= 7 and
            len(box_clean)/len(box_total) >= 0.8 and
            len(box_cd73_dim) > 0 and len(box_cd73_bright) > 0
          ):

            if annotate_points:

                # stats_dirate current box with grid coordinates
                label = '*'  # (e, j)

                for ax in [ax_adeno_img]:

                    ax.annotate(
                        label, size=0.8,
                        xy=(col_start, row_start),
                        xytext=(1.0, -2.0),
                        textcoords='offset points',
                        ha='left', va='bottom', weight='normal',
                        alpha=1.0, color='lime',
                        bbox=dict(boxstyle='round,pad=0.1',
                                  fc='yellow', alpha=0.0)
                        )

            # calculate metabolite signals
            if density:
                # calculate per cell metabolite concentrations
                concentration_adeno = (
                    np.mean(box_adeno)/len(box_total)
                    )
            else:
                # calculate mean metabolite concentration
                concentration_adeno = np.mean(box_adeno)

            # store current measurements
            measurements = [
                [e, j, float(len(box_cd73_bright)),
                 float(len(box_cd73_dim)),
                 float(concentration_adeno)]
                ]
            box_measurements = pd.DataFrame(
                measurements, columns=['row', 'column',
                                       'cd73_bright',
                                       'cd73_dim',
                                       'concentration_adeno']
                )

            # add box measurements overall measurements dataframe
            overall_measurements = pd.concat(
                [overall_measurements, box_measurements],
                ignore_index=True
                )

        col_start = col_end

    row_start = row_end

# add cell centroids to figure
for ax in [ax_adeno_img]:

    ax.scatter(
        cd73_dim['X_centroid'], cd73_dim['Y_centroid'],
        s=0.4, c='tab:blue', lw=0.0, zorder=1, marker='o',
        alpha=0.2
        )

    ax.scatter(
        cd73_bright['X_centroid'], cd73_bright['Y_centroid'],
        s=0.4, c='r', ec='k', lw=0.0, zorder=1, marker='o',
        alpha=0.2
        )

# add plots to figure
for name, ax, value in zip(
  ['adenosine'],
  [ax_adeno_plot],
  ['concentration_adeno']
  ):

    overall_measurements['ratio'] = (
        overall_measurements['cd73_bright'] /
        overall_measurements['cd73_dim']
        )

    # categorize boxes according to their ratio of
    # CD73 bright to CD73 dim cells
    overall_measurements['category'] = [
        '>0.9' if i > 0.9 else '0.9-0.3' if
        0.3 <= i <= 0.9 else '<0.3' for
        i in overall_measurements['ratio']
        ]

    # make corresponding numerical category codes for plotting
    overall_measurements['annotate'] = [
        0.0 if i == '>0.9' else 1.0 if i == '0.9-0.3'
        else 2.0 for i in overall_measurements['category']
        ]

    # add stripplots to figure
    sns.stripplot(
        x='category', y=value,
        data=overall_measurements, color='k', s=0.75,
        edgecolor='w', linewidth=0.05, ax=ax
        )

    # plot mean line
    sns.boxplot(showmeans=True,
                meanline=True,
                meanprops={'color': 'k', 'ls': '-', 'lw': 0.5},
                medianprops={'visible': False},
                whiskerprops={'visible': False},
                zorder=10,
                x='category',
                y=value,
                data=overall_measurements,
                showfliers=False,
                showbox=False,
                showcaps=False,
                ax=ax)

    ax.set_xlabel(
        '$\mathregular{CD73^{bright}/CD73^{dim}}$ cells per box',
        fontsize=6, labelpad=8.0
        )

    if density:
        y_label = f'Density-corrected [{name}]'
    else:
        y_label = f'mean {name} intensity per box'

    ax.set_ylabel(y_label, fontsize=6, labelpad=10.0)
    ax.tick_params(axis='both', which='major', labelsize=6)

    # extract boxes containing a high ratio of bright-to-dim cells
    cd73_high_boxes = overall_measurements[
        overall_measurements['category'] == '>0.9'
        ]

    # extract boxes containing an intermediate ratio of
    # bright-to-dim cells
    cd73_medium_boxes = overall_measurements[
        overall_measurements['category'] == '0.9-0.3'
        ]

    # extract boxes containing a low ratio of bright-to-dim cells
    cd73_low_boxes = overall_measurements[
        overall_measurements['category'] == '<0.3'
        ]

    # identify ratio category associated with fewest boxes
    # to get n for random sampling
    sample_size = min(
        [len(cd73_high_boxes),
         len(cd73_medium_boxes),
         len(cd73_low_boxes)]
         )

    # take equal number of boxes from each ratio category
    high_sample = (
        cd73_high_boxes[value]
        .sample(n=sample_size, random_state=3)
        .reset_index(drop=True)
        .rename('>0.9')
        )
    medium_sample = (
        cd73_medium_boxes[value]
        .sample(n=sample_size, random_state=3)
        .reset_index(drop=True)
        .rename('0.9-0.3')
        )
    low_sample = (
        cd73_low_boxes[value]
        .sample(n=sample_size, random_state=3)
        .reset_index(drop=True)
        .rename('<0.3')
        )

    samples = pd.concat(
        [high_sample, medium_sample, low_sample], axis=1
        )

    # compute one-way ANOVA
    fval, pval = f_oneway(
        samples['<0.3'],
        samples['0.9-0.3'],
        samples['>0.9'],
        )

    # perform multiple pairwise comparisons (Tukey's HSD)
    samples_melt = pd.melt(
        samples.reset_index(), id_vars=['index'],
        value_vars=['>0.9', '0.9-0.3', '<0.3']
        )

    samples_melt.columns = ['index', 'treatments', 'value']

    res = stat()
    res.tukey_hsd(
        df=samples_melt, res_var='value',
        xfac_var='treatments',
        anova_model='value ~ C(treatments)')

    # add statistics to figure
    ax.text(
        x=0.48, y=0.94,
        s=f"F={round(fval, 2)}, " +
        f"pval={round(pval, 4)}",
        horizontalalignment='left', verticalalignment='center',
        size='xx-small', transform=ax.transAxes,
        color='black', weight='normal')

# title figure
fig.suptitle('brain6', x=0.20, y=0.94, fontsize=9)

# save figure
plt.savefig(os.path.join(save_dir, 'brain6.png'), dpi=1000)
plt.close('all')

# save Tukey's HSD summary table
res.tukey_summary.to_csv(
    os.path.join(save_dir, 'tukey_hsd.csv'), index=False
    )
