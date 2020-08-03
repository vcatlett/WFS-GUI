import sys
import os
import time
from subprocess import Popen, PIPE
import matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from astropy.io import fits
from PyQt5.QtCore import Qt, QTimer
from scipy.stats import multivariate_normal
from PyQt5.QtWidgets import (QWidget, QPushButton, QMessageBox, QVBoxLayout, QGridLayout, QApplication,
                             QMainWindow, QCheckBox, QDial, QTabWidget, QSlider, QRadioButton, QLabel, QComboBox,
                             QTableWidget, QTableWidgetItem, QHeaderView)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
matplotlib.use('Qt5Agg')

os.environ['LD_LIBRARY_PATH'] = '/usr/local/lib:/usr/lib:/usr/local/lib64:/usr/lib64'

print('\n Loading... \n')

class importDefaults:
    def importConfig(self):
        global plateScales, scopeRadius, ps

        plateScales = [1.515E-6, 1.515E-6, 1.515E-6, 1.515E-6, 3.030E-6]  # 1, 2, 4, 8, 16
        scopeRadius = 1.1E9  # nm
        ps = plateScales[4]

    def importMatrices(self):
        global z20_xs, zcoeffs, lis, z16, z8, z4, z2, z1, rotMatrix_16, rotMatrix_8, rotMatrix_4, rotMatrix_2
        z20_xs = np.linspace(1, 20, 20)
        zcoeffs = list(np.zeros(20))
        lis = [0, 0.5, 0.9, 0.99, 0.999, 1.0]

        #zpath = '/Users/victoriacatlett/Desktop/2020_Summer/REU/Matrices/Zernike/'
        #rpath = '/Users/victoriacatlett/Desktop/2020_Summer/REU/Matrices/Rotation/'
        zpath = './matrices/Zernike/'
        rpath = './matrices/Rotation/'

        z16 = pd.read_table(zpath + 'Z20_16.txt', header=None, names=None).values
        z8 = pd.read_table(zpath + 'Z20_8.txt', header=None, names=None).values
        z4 = pd.read_table(zpath + 'Z14_4.txt', header=None, names=None).values
        z2 = pd.read_table(zpath + 'Z5_2.txt', header=None, names=None).values
        z1 = np.zeros((1, 2))
        rotMatrix_16 = np.transpose(pd.read_table(rpath + 'Rotation_Reconstruction_16.txt', names=None).values)
        rotMatrix_8 = np.transpose(pd.read_table(rpath + 'Rotation_Reconstruction_8.txt', names=None).values)
        rotMatrix_4 = np.transpose(pd.read_table(rpath + 'Rotation_Reconstruction_4.txt', names=None).values)
        rotMatrix_2 = np.transpose(pd.read_table(rpath + 'Rotation_Reconstruction_2.txt', names=None).values)

    def importPlotDefaults(self):
        global gla, gca, fca, sl, lowerSlide, upperSlide, image_mode, boxsize, back_color, image_data, dark_data, dCenters
        gla = 0
        gca = 0
        fca = 0
        sl = 0
        lowerSlide = -0.25
        upperSlide = 2.0
        image_mode = 16  # 8, 4, 2, 1
        boxsize = 8
        back_color = [0.195, 0.195, 0.195]
        dCenters = np.linspace(3.5, 123.5, 16)
        image_data = np.zeros((128, 128))
        dark_data = np.random.normal(loc=0.0, scale=1.0, size=(128, 128)) * 0.05

    def dataDefaults(self):
        global centMu, centSigma, dotCovar, dotAmp, rms, theta, li_coeff
        centMu = 3.5
        centSigma = 0.5
        dotCovar = 1.0
        dotAmp = 10
        rms = 0.0
        theta = 0.0
        li_coeff = 0.0

    def importBools(self):
        global noDark, dark_sub, taking_dark
        noDark = True
        dark_sub = False
        taking_dark = False


class NUVU:
    def getNPS(self):
        stream = os.popen('nps nps_1')
        output = stream.readlines()
        tbl = []
        for i in range(10, 26, 1):
            line = output[i].split()
            line = [x for x in line if (x != '|' and x != 'sec')]
            tbl.append(line)

        return tbl

    def connect(self):
        global p
        print('\n Attempting to connect... \n')
        p = Popen('./NUVU/bin/wfs_test', shell=True, stdout=PIPE, stdin=PIPE)
        value = 'a' + '\n'
        value = bytes(value, 'UTF-8')
        p.stdin.write(value)
        '''print('Printing errors...')
        for i in range(5):
            result = p.stdout.readline().strip()
            print(result.decode('ascii'))'''

    def setNames(self):
        value = 'n' + '\n'
        value = bytes(value, 'UTF-8')
        p.stdin.write(value)

        value = 'img' + '\n'
        value = bytes(value, 'UTF-8')
        p.stdin.write(value)

        value = 'dark' + '\n'
        value = bytes(value, 'UTF-8')
        p.stdin.write(value)

        print('\n Default image names set \n')

    def openShutter(self):
        print('\n Attempting to open shutter... \n')
        value = 'o' + '\n'
        value = bytes(value, 'UTF-8')
        p.stdin.write(value)

    def closeShutter(self):
        print('\n Attempting to close shutter... \n')
        value = 's' + '\n'
        value = bytes(value, 'UTF-8')
        p.stdin.write(value)

    def singleCapture(self):
        #print('\n Taking single image... \n')
        try:
            value = 'b' + '\n'
            value = bytes(value, 'UTF-8')
            p.stdin.write(value)

            '''for i in range(1):
                result = p.stdout.readline().strip()
                print(result.decode('ascii'))'''

            #value = 'img' + '\n'
            #value = bytes(value, 'UTF-8')
            #p.stdin.write(value)

            #value = '1' + '\n'
            #value = bytes(value, 'UTF-8')
            #p.stdin.write(value)

        except (BrokenPipeError, IOError):
            print('\n PIPE BROKE \n')

    def contAcq(self):
        print('\n Continuous capture... \n')
        value = 'e' + '\n'
        value = bytes(value, 'UTF-8')
        p.stdin.write(value)
        #p.stdin.flush()

        value = '1' + '\n'
        value = bytes(value, 'UTF-8')
        p.stdin.write(value)
        #p.stdin.flush()

        value = 's' + '\n'
        value = bytes(value, 'UTF-8')
        p.stdin.write(value)
        #p.stdin.flush()

        for i in range(243):
            result = p.stdout.readline().strip()
            print(result.decode('ascii'))
        print('\n Continuous capture successful \n')

    def darkCapture(self):
        #print('\n Taking Dark Frame... \n')
        try:
            value = 'j' + '\n'
            value = bytes(value, 'UTF-8')
            p.stdin.write(value)

            value = 'dark' + '\n'
            value = bytes(value, 'UTF-8')
            p.stdin.write(value)

            value = '1' + '\n'
            value = bytes(value, 'UTF-8')
            p.stdin.write(value)

        except (BrokenPipeError, IOError):
            print('\n PIPE BROKE \n')

        #print('Finished dark commands')

    def disconnect(self):
        print('\n Disconnecting from WFS... \n')
        value = 'h' + '\n'
        value = bytes(value, 'UTF-8')
        p.stdin.write(value)
        #p.stdin.flush()

        #for i in range(2):
        #result = p.stdout.readlines()
        #print(result.decode('ascii'))

    def close(self):
        print('\n Closing connection program... \n')
        value = 'x' + '\n'
        value = bytes(value, 'UTF-8')
        p.stdin.write(value)
        #p.stdin.flush()
        #for i in range(2):
            #result = p.stdout.readline().strip()
            #print(result.decode('ascii'))
        print('\n Closed \n')


class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, width, height, dpi, color, size, pad, parent=None):
        fig = plt.figure(figsize=(width, height), dpi=dpi, facecolor=color)
        self.axes = fig.add_axes((1-size, 1-size, size-pad, size-pad), frameon=False)
        self.axes.axes.tick_params(colors=back_color)
        super(MplCanvas, self).__init__(fig)


class UI_Popups:
    def noDark(self):
        msg = QMessageBox()
        msg.setWindowTitle('Dark Warning')
        msg.setText('You have not taken a dark yet.')
        msg.exec_()

    def oldDark(self):
        msg = QMessageBox()
        msg.setWindowTitle('Dark Warning')
        msg.setText('Your last dark was too long ago.')


class makeData:
    def newData(self):
        print('\n Creating New Image... \n')
        # data = np.zeros((128, 128))
        boxsize = int(128/image_mode)
        boxCenters = np.linspace(centMu, 127-centMu, image_mode)
        boxcentxs, boxcentys = np.meshgrid(boxCenters, boxCenters)

        data = np.random.normal(loc=0.0, scale=1.0, size=(128, 128))*0.01
        centxs = np.random.normal(centMu, centSigma, size=(16, 16))
        centys = np.random.normal(centMu, centSigma, size=(16, 16))
        shift = 63.5

        theta = np.random.normal(loc=0, scale=1)*0.75
        c = np.cos(np.deg2rad(theta))
        s = np.sin(np.deg2rad(theta))
        rot_by_theta = np.array(((c, -s), (s, c)))

        for i in range(image_mode):
            indx0i = boxsize * i
            indx1i = indx0i + boxsize
            for j in range(image_mode):
                indx0j = boxsize * j
                indx1j = indx0j + boxsize

                box = np.zeros((boxsize, boxsize), dtype=np.float64)
                centx = centxs[i][j]
                centy = centys[i][j]
                pre_rot = [centx - centMu + boxcentxs[i][j] - shift, centy - centMu + boxcentys[i][j] - shift]
                rot_xy = np.matmul(rot_by_theta, pre_rot)
                rot_xy[0] = rot_xy[0] - boxcentxs[i][j] + shift + centMu
                rot_xy[1] = rot_xy[1] - boxcentys[i][j] + shift + centMu
                var = multivariate_normal(mean=rot_xy, cov=[[dotCovar*(16/image_mode)**2, 0], [0, dotCovar*(16/image_mode)**2]])
                rad = np.sqrt((boxcentxs[i][j] - shift) ** 2 + (boxcentys[i][j] - shift) ** 2)

                for k in range(boxsize):
                    for l in range(boxsize):
                        if rad < 69:
                            boxfill = var.pdf([k, l]) * dotAmp * (16/image_mode)**2
                            box[l][k] = boxfill
                        else:
                            box[l][k] = 0

                data[indx0i:indx1i, indx0j:indx1j] += box
        return data


class figCalcs:

    def get_centroids(self, img):
        global dCenters, rms, theta
        boxsize = int(128/image_mode)
        dCenters = np.linspace(centMu, 127-centMu, image_mode)

        fcentxs = np.zeros((image_mode, image_mode))
        fcentys = np.zeros((image_mode, image_mode))
        xfromcent = np.zeros((image_mode, image_mode))
        yfromcent = np.zeros((image_mode, image_mode))
        disps = np.zeros((image_mode, image_mode))

        for i in range(image_mode):
            indx0i = boxsize * i
            for j in range(image_mode):
                rad = np.sqrt((dCenters[i]-63.5)**2 + (dCenters[j]-63.5)**2)
                indx0j = boxsize * j
                box_ij = self.pull_box(img, indx0i, indx0j, boxsize)

                if rad < 69:
                    centroid = self.box_centroid(boxsize, box_ij)
                else:
                    centroid = [0, 0, 0]

                fcentxs[i][j] = centroid[0] + centMu
                fcentys[i][j] = centroid[1] + centMu
                xfromcent[i][j] = centroid[0]
                yfromcent[i][j] = centroid[1]
                disps[i][j] = centroid[2]

        flat_xy = np.array(list(xfromcent.flatten('C')) + list(yfromcent.flatten('C')))
        z_coeffs = self.calcZernike(flat_xy)
        theta = self.calcRotation(flat_xy)
        rms = self.calcRMS(disps.flatten('C'))

        return fcentxs, fcentys, disps, z_coeffs

    def pull_box(self, data, iIndx, jIndx, boxsize):
        box = np.zeros((boxsize, boxsize), dtype=np.float64)
        for k in range(boxsize):
            for l in range(boxsize):
                box[k][l] = data[iIndx + k][jIndx + l]
        return box

    def box_centroid(self, boxsize, boxdata):
        x_num = 0
        x_den = 0
        y_num = 0
        y_den = 0

        for i in range(int(boxsize)):
            x_num = x_num + (i - centMu) * sum(boxdata[:, i])
            x_den = x_den + sum(boxdata[:, i])
            y_num = y_num + (i - centMu) * sum(boxdata[i, :])
            y_den = y_den + sum(boxdata[i, :])

        if x_den > 0 and y_den > 0:
            cent_x = x_num / x_den
            cent_y = y_num / y_den
        else:
            cent_x = 0
            cent_y = 0

        disp = np.sqrt(cent_x ** 2 + cent_y ** 2)
        return cent_x, cent_y, disp

    def calcRMS(self, disps):
        n = len(disps)
        sqrs = [d ** 2 for d in disps]
        rms = np.sqrt(sum(sqrs) / n)
        return rms

    def calcZernike(self, flat_xy):
        if image_mode == 16:
            z_coeffs = [z*ps*scopeRadius for z in np.matmul(z16, flat_xy)]
        elif image_mode == 8:
            z_coeffs = [z*ps*scopeRadius for z in np.matmul(z8, flat_xy)]
        elif image_mode == 4:
            z_coeffs = [z*ps*scopeRadius for z in np.matmul(z4, flat_xy)]
            z_coeffs.extend(np.zeros(6))
        elif image_mode == 2:
            z_coeffs = [z*ps*scopeRadius for z in np.matmul(z2, flat_xy)]
            z_coeffs.extend(np.zeros(15))
        #elif image_mode == 1:
           # z_coeffs = list(flat_xy)
           # z_coeffs.extend(np.zeros(18))
        else:
            z_coeffs = list(np.zeros(20))
        return z_coeffs

    def calcRotation(self, flat_xy):
        if image_mode == 16:
            rot_angle = np.matmul(rotMatrix_16, flat_xy)
        elif image_mode == 8:
            rot_angle = np.matmul(rotMatrix_8, flat_xy)
        elif image_mode == 4:
            rot_angle = np.matmul(rotMatrix_4, flat_xy)
        elif image_mode == 2:
            rot_angle = np.matmul(rotMatrix_2, flat_xy)
        else:
            rot_angle = 0
        return rot_angle


class saveImage:
    def savePNG(self, fig):
        #plot = fig.axes
        sfig = fig.figure
        plt.margins(0, 0)
        #sfig.savefig('/Users/victoriacatlett/Desktop/test.png', aspect=1.0)
        timestr = time.strftime("_%Y_%m_%d-%H_%M_%S")
        sfig.savefig('./PNGs/img'+timestr+'.png', aspect=1.0)
        print('\n Saved PNG Image img'+timestr+'.png \n')

    def saveFITS(self):
        timestr = time.strftime("_%Y_%m_%d-%H_%M_%S")
        hdu = fits.PrimaryHDU(data=image_data)
        hdu.writeto('./FITS/img'+timestr+'.fits')
        print('\n Saved FITS Image img'+timestr+'.png \n')


class table():
    def makeTable(self):
        data = NUVU().getNPS()
        shape = np.shape(data)
        tbl = QTableWidget(shape[0], shape[1]-1, parent=None)
        tbl.setHorizontalHeaderLabels(['Name', 'State', 'Boot Delay (sec)', 'Default'])
        hdr = tbl.horizontalHeader()
        hdr.setSectionResizeMode(0, QHeaderView.Stretch)
        hdr.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        hdr.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        hdr.setSectionResizeMode(3, QHeaderView.ResizeToContents)
        for i in range(shape[0]):
            for j in range(1, shape[1], 1):
                item = QTableWidgetItem()
                item.setText(data[i][j])
                tbl.setItem(i, j-1, item)

        return tbl


class MainWindow(QMainWindow, QApplication):

    def __init__(self, *args, **kwargs):
        # These lines create the window and import/load everything
        super(MainWindow, self).__init__(*args, **kwargs)
        print('\n Gathering data... \n')
        self.importAll()
        print('\n Connecting to WFS... \n')
        NUVU().connect()
        print('\n Opening shutter... \n')
        NUVU().openShutter()
        print('\n Setting default image names... \n')
        NUVU().setNames()
        print('\n Creating window... \n')
        win = self.makeWindow()
        win.show()

        # These run when the window closes
        app.aboutToQuit.connect(self.closeEvent)
        sys.exit(app.exec_())

    def closeEvent(self):
        print('\n Exiting... \n')
        NUVU().disconnect()
        print('Disconnected')

    def importAll(self):
        importDefaults().importConfig()
        importDefaults().importMatrices()
        importDefaults().importPlotDefaults()
        importDefaults().dataDefaults()
        importDefaults().importBools()

    def makeWindow(self):
        win = QWidget()
        #win.setStyleSheet("background-color:rgb(50, 50, 50);")

        sc1 = MplCanvas(5.0, 4.0, 100, back_color, 1, 0)
        sc2 = MplCanvas(4.0, 2.0, 100, back_color, 0.85, 0.05)

        box_gl, box_gc, box_fc, box_dark, refresh_btn, dark_btn, background_btn, png_btn, fits_btn, ss_btn \
            = self.makeImageButtons()

        boxSlider = self.makeSlider(0, 128)
        liSlider = self.makeSlider(0, len(lis)-1)
        rot_dial = self.makeDial()
        rms_lbl = self.makeRMSLabel()
        rot_lbl = self.makeRotLabel()
        box_lbl = QLabel('%i px' % sl)
        li_lbl = QLabel(str(li_coeff))

        modeCB = QComboBox()
        modeCB.addItems(['1x1', '2x2', '4x4', '8x8', '16x16'])
        modeCB.setCurrentIndex(4)
        self.connectButtons(box_gl, box_gc, box_fc, box_dark, refresh_btn, dark_btn, png_btn, fits_btn, ss_btn,
                            modeCB, boxSlider, liSlider, sc1, sc2, rot_dial, rms_lbl, rot_lbl, box_lbl, li_lbl)

        tab1_layout = self.makeTab1(refresh_btn, dark_btn, background_btn, box_dark)
        tab2_layout = self.makeTab2(boxSlider, box_lbl, liSlider, li_lbl, modeCB)
        tab3_layout = self.makeTab3(png_btn, fits_btn, ss_btn)
        tab4_layout = self.makeTab4()

        tabs = self.makeAllTabs(tab1_layout, tab2_layout, tab3_layout, tab4_layout)

        lbox = self.makeLbox(tabs, sc2, rot_dial, box_gc, box_fc, box_gl, rms_lbl, rot_lbl)
        rbox = self.makeRbox(sc1)

        hbox = QGridLayout()
        hbox.setColumnStretch(0, 4.0)
        hbox.setColumnStretch(1, 6.0)

        hbox.addLayout(lbox, 0, 0)
        hbox.addLayout(rbox, 0, 1)

        win.setLayout(hbox)
        win.setWindowTitle("WFS Alignment GUI")
        win.setGeometry(100, 100, 1200, 700)
        win.setStyleSheet(open('./styles.css').read())

        return win

    '''def makeModeButtons(self):
        rad16 = QRadioButton("16 x 16")
        rad8 = QRadioButton("8 x 8")
        rad4 = QRadioButton("4 x 4")
        rad2 = QRadioButton("2 x 2")
        rad1 = QRadioButton("1 x 1")
        rad16.setChecked(True)

        return rad1, rad2, rad4, rad8, rad16'''

    def makeImageButtons(self):
        box_gl = QCheckBox("Grid Lines")
        box_gl.setChecked(False)
        box_gl.setStyleSheet(open('./styles.css').read())

        box_gc = QCheckBox("Grid Centers")
        box_gc.setChecked(False)
        box_gc.setStyleSheet(open('./styles.css').read())

        box_fc = QCheckBox("Found Centers")
        box_fc.setChecked(False)
        box_fc.setStyleSheet(open('./styles.css').read())

        box_dark = QCheckBox("Dark Subtraction")
        box_dark.setChecked(False)
        box_dark.setStyleSheet(open('./styles.css').read())

        refresh_btn = QPushButton('Refresh Image')
        dark_btn = QPushButton('Take Dark')
        background_btn = QPushButton('Take Background')
        png_btn = QPushButton('Save PNG')
        fits_btn = QPushButton('Save FITS')
        ss_btn = QPushButton('Take Screenshot')

        return box_gl, box_gc, box_fc, box_dark, refresh_btn, dark_btn, background_btn, png_btn, fits_btn, ss_btn

    def makeSlider(self, sMin, sMax):
        slider = QSlider(Qt.Horizontal)
        slider.setMinimum(sMin)
        slider.setMaximum(sMax)
        slider.setStyleSheet(open('./styles.css').read())
        #slider.setTickPosition(Qt.TicksBelow)

        return slider

    def makeDial(self):
        rot_dial = QDial()
        rot_dial.setRange(0, 360)  # 0: -3, 60: -2, 120: -1, ...
        rot_dial.setValue(180)
        rot_dial.setStyleSheet(open('./styles.css').read())

        return rot_dial

    def makeRMSLabel(self):
        rms_lbl = QLabel()
        rms_lbl.setText('RMS = None')
        rms_lbl.setStyleSheet(open('./styles.css').read())

        return rms_lbl

    def makeRotLabel(self):
        rot_lbl = QLabel()
        rot_lbl.setText('Theta = None')
        rot_lbl.setStyleSheet(open('./styles.css').read())
        return rot_lbl

    def connectButtons(self, box_gl, box_gc, box_fc, box_dark, refresh_btn, dark_btn, png_btn, fits_btn, ss_btn, modeCB,
                       boxSlider, liSlider, sc1, sc2, rot_dial, rms_lbl, rot_lbl, box_lbl, li_lbl):
        boxSlider.valueChanged.connect(lambda: self.boxSliderState(boxSlider, box_lbl, rot_dial, sc1, sc2))
        liSlider.sliderReleased.connect(lambda: self.liSliderState(liSlider, li_lbl, rot_dial, sc1, sc2))
        refresh_btn.clicked.connect(lambda: self.refreshData(rot_dial, sc1, sc2, rms_lbl, rot_lbl))
        dark_btn.clicked.connect(lambda: self.takeDark(rot_dial, sc1, sc2))
        png_btn.clicked.connect(lambda: saveImage.savePNG(self, sc1))
        fits_btn.clicked.connect(lambda: saveImage.saveFITS(self))

        box_gl.stateChanged.connect(lambda: self.boxstate(box_gl, box_gc, box_fc, box_dark, sc1, sc2, rot_dial))
        box_gc.stateChanged.connect(lambda: self.boxstate(box_gl, box_gc, box_fc, box_dark, sc1, sc2, rot_dial))
        box_fc.stateChanged.connect(lambda: self.boxstate(box_gl, box_gc, box_fc, box_dark, sc1, sc2, rot_dial))
        box_dark.stateChanged.connect(lambda: self.boxstate(box_gl, box_gc, box_fc, box_dark, sc1, sc2, rot_dial))
        modeCB.currentIndexChanged.connect(lambda: self.modeChange(modeCB, sc1, sc2, rot_dial, rms_lbl, rot_lbl))

    def makeTab1(self, refresh_btn, dark_btn, background_btn, box_dark):
        tab1_layout = QGridLayout()
        tab1_layout.addWidget(refresh_btn, 0, 0)
        tab1_layout.addWidget(dark_btn, 0, 1)
        tab1_layout.addWidget(background_btn, 0, 2)
        tab1_layout.addWidget(box_dark, 1, 1)

        '''tab1_layout.addWidget(QLabel(text='1x1 Box Side Length'), 1, 0)
        tab1_layout.addWidget(boxSlider, 1, 1)
        tab1_layout.addWidget(box_lbl, 1, 2)

        tab1_layout.addWidget(QLabel(text='LI Coefficient'), 2, 0)
        tab1_layout.addWidget(liSlider, 2, 1)
        tab1_layout.addWidget(li_lbl, 2, 2)'''

        return tab1_layout

    def makeTab2(self, boxSlider, box_lbl, liSlider, li_lbl, modeCB):
        tab2_layout = QGridLayout()

        tab2_layout.addWidget(QLabel(text='Lenslet Mode'), 0, 0)
        tab2_layout.addWidget(modeCB, 0, 1)

        tab2_layout.addWidget(QLabel(text='1x1 Box Side Length'), 1, 0)
        tab2_layout.addWidget(boxSlider, 1, 1)
        tab2_layout.addWidget(box_lbl, 1, 2)

        tab2_layout.addWidget(QLabel(text='LI Coefficient'), 2, 0)
        tab2_layout.addWidget(liSlider, 2, 1)
        tab2_layout.addWidget(li_lbl, 2, 2)

        return tab2_layout

    def makeTab3(self, png_btn, fits_btn, ss_btn):
        tab3_layout = QGridLayout()
        tab3_layout.addWidget(png_btn, 1, 0)
        tab3_layout.addWidget(fits_btn, 1, 1)
        tab3_layout.addWidget(ss_btn, 1, 2)
        return tab3_layout

    def makeTab4(self):
        tab4_layout = QGridLayout()

        '''data = NUVU().getNPS()
        shape = np.shape(data)
        tbl = QTableWidget(shape[0], shape[1], parent=None)
        tbl.setHorizontalHeaderLabels(['Plug', 'Name', 'State', 'Boot Delay (sec)', 'Default'])

        for i in range(shape[0]):
            for j in range(shape[1]):
                tbl.setItem(i, j, data[i][j])'''
        tbl = table().makeTable()
        tab4_layout.addWidget(tbl)

        return tab4_layout

    def makeAllTabs(self, tab1Layout, tab2Layout, tab3Layout, tab4Layout):
        tabs = QTabWidget()
        tabs.setStyleSheet(open('./styles.css').read())
        tab1 = QWidget()
        tab2 = QWidget()
        tab3 = QWidget()
        tab4 = QWidget()
        tab1.setLayout(tab1Layout)
        tab2.setLayout(tab2Layout)
        tab3.setLayout(tab3Layout)
        tab4.setLayout(tab4Layout)

        tabs.addTab(tab1, "Image Controls")
        tabs.addTab(tab2, "Camera Controls")
        tabs.addTab(tab3, "Save Images")
        tabs.addTab(tab4, "NPS Info")

        return tabs

    def makeLbox(self, tabs, sc2, rot_dial, box_gc, box_fc, box_gl, rms_lbl, rot_lbl):
        lbox = QGridLayout()

        lbox.setRowStretch(0, 3)
        lbox.setRowStretch(1, 4)
        lbox.setRowStretch(2, 3)
        #lbox.setRowStretch(3, 1)

        lbox_0 = QGridLayout()
        lbox_1 = QGridLayout()
        lbox_1_0 = QGridLayout()
        lbox_1_1 = QGridLayout()
        lbox_1_2 = QGridLayout()
        lbox_2 = QGridLayout()
        lbox_3 = QGridLayout()

        #lbox_2.setColumnStretch(0, 1)
        #lbox_2.setColumnStretch(1, 9)

        lbox_1_0.addWidget(box_gl, 0, 0)
        lbox_1_0.addWidget(box_gc, 1, 0)
        lbox_1_0.addWidget(box_fc, 2, 0)

        lbox_1_1.addWidget(rot_dial, 0, 0)

        lbox_1_2.addWidget(rms_lbl, 0, 0)
        lbox_1_2.addWidget(rot_lbl, 1, 0)

        #x_placeholder = QLabel(text='X Placeholder')
        #y_placeholder = QLabel(text='Y Placeholder')

        lbox_0.addWidget(tabs)
        lbox_1.addLayout(lbox_1_0, 0, 0)
        lbox_1.addLayout(lbox_1_1, 0, 1)
        lbox_1.addLayout(lbox_1_2, 0, 2)
        #lbox_2.addWidget(y_placeholder, 0, 0)
        lbox_2.addWidget(sc2, 0, 0)
        #lbox_3.addWidget(x_placeholder)

        lbox.addLayout(lbox_0, 0, 0)
        lbox.addLayout(lbox_1, 1, 0)
        lbox.addLayout(lbox_2, 2, 0)
        lbox.addLayout(lbox_3, 3, 0)

        return lbox

    def makeRbox(self, sc1):
        rbox = QGridLayout()
        rbox.addWidget(sc1)
        return rbox

    def boxstate(self, b1, b2, b3, b4, sc1, sc2, rot_dial):
        global gla, gca, fca, dark_sub

        if b1.isChecked():
            gla = 1
        else:
            gla = 0

        if b2.isChecked():
            gca = 1
        else:
            gca = 0

        if b3.isChecked():
            fca = 1
        else:
            fca = 0

        if b4.isChecked():
            if noDark:
                UI_Popups().noDark()
                b4.setChecked(False)
                dark_sub = False
            else:
                dark_sub = True
        else:
            dark_sub = False

        self.update_plot(rot_dial, sc1, sc2)

    def radiostate(self, r, sc1, sc2, rot_dial, rms_lbl, rot_lbl):
        global image_mode, centMu, centSigma, boxsize
        boxsize = int(128 / image_mode)

        if r.text() == "16 x 16":
            if r.isChecked() == True:
                image_mode = 16
                centMu = 3.5
                centSigma = 0.5
                self.refreshData(rot_dial, sc1, sc2, rms_lbl, rot_lbl)

        if r.text() == "8 x 8":
            if r.isChecked():
                image_mode = 8
                centMu = 7.5
                centSigma = 1.0
                self.refreshData(rot_dial, sc1, sc2, rms_lbl, rot_lbl)

        if r.text() == "4 x 4":
            if r.isChecked():
                image_mode = 4
                centMu = 15.5
                centSigma = 2.0
                self.refreshData(rot_dial, sc1, sc2, rms_lbl, rot_lbl)

        if r.text() == "2 x 2":
            if r.isChecked():
                image_mode = 2
                centMu = 31.5
                centSigma = 4.0
                self.refreshData(rot_dial, sc1, sc2, rms_lbl, rot_lbl)

        if r.text() == "1 x 1":
            if r.isChecked():
                image_mode = 1
                centMu = 63.5
                centSigma = 8.0
                self.refreshData(rot_dial, sc1, sc2, rms_lbl, rot_lbl)

    def modeChange(self, cb, sc1, sc2, rot_dial, rms_lbl, rot_lbl):
        global image_mode, centMu, centSigma, boxsize, zcoeffs

        # TO DO: Change to dictionary
        boxsize = int(128 / image_mode)
        zcoeffs = list(np.zeros(20))

        if cb.currentIndex() == 4:
            image_mode = 16
            centMu = 3.5
            centSigma = 0.5
            self.refreshData(rot_dial, sc1, sc2, rms_lbl, rot_lbl)

        if cb.currentIndex() == 3:
            image_mode = 8
            centMu = 7.5
            centSigma = 1.0
            self.refreshData(rot_dial, sc1, sc2, rms_lbl, rot_lbl)

        if cb.currentIndex() == 2:
            image_mode = 4
            centMu = 15.5
            centSigma = 2.0
            self.refreshData(rot_dial, sc1, sc2, rms_lbl, rot_lbl)

        if cb.currentIndex() == 1:
            image_mode = 2
            centMu = 31.5
            centSigma = 4.0
            self.refreshData(rot_dial, sc1, sc2, rms_lbl, rot_lbl)

        if cb.currentIndex() == 0:
            image_mode = 1
            centMu = 63.5
            centSigma = 8.0
            self.refreshData(rot_dial, sc1, sc2, rms_lbl, rot_lbl)

    def boxSliderState(self, boxSlider, box_lbl, rot_dial, sc1, sc2):
        global sl
        sl = boxSlider.value()
        box_lbl.setText('%i px' % sl)
        self.update_plot(rot_dial, sc1, sc2)
        return boxSlider

    def liSliderState(self, liSlider, li_lbl, rot_dial, sc1, sc2):
        global li_coeff
        li_coeff = lis[liSlider.value()]
        li_lbl.setText(str(li_coeff))
        self.update_plot(rot_dial, sc1, sc2)
        return liSlider

    def refreshData(self, rot_dial, sc1, sc2, rms_lbl, rot_lbl):
        global image_data, dark_data, taking_dark
        taking_dark = False
        print('\n Taking frame...')
        #image_data = makeData().newData()
        avg_data = np.zeros((130, 128))
        #NUVU().openShutter()
        for i in range(100):
            NUVU().singleCapture()
            newImg = fits.getdata('./NUVU/bin/img.fits')
            avg_data = np.add(avg_data, newImg)
            print('i ' + str(newImg[7][7]))
        print('\n')
        image_data = avg_data/100
        print('Overall: ' + str(image_data[7][7]))
        image_data = np.delete(image_data, 1, 0)
        image_data = np.delete(image_data, 0, 0)
        self.update_plot(rot_dial, sc1, sc2)
        rms_lbl.setText('RMS = %0.2f px' % rms)
        rot_lbl.setText('Theta = %0.2f degrees' % theta)

    def timer(self, rot_dial, sc1, sc2, rms_lbl, rot_lbl):
        self.timer = QTimer()
        self.timer.setInterval(100)
        self.timer.timeout.connect(lambda: self.refreshData(rot_dial, sc1, sc2, rms_lbl, rot_lbl))
        self.timer.start()

    def takeDark(self, rot_dial, sc1, sc2):
        global image_data, dark_data, taking_dark, noDark
        print('\n Taking Dark Frames...')
        taking_dark = True
        noDark = False
        dark_data = np.zeros((130, 128))
        NUVU().closeShutter()
        for i in range(100):
            #print('Dark ' + str(i+1))
            NUVU().darkCapture()
            newDark = fits.getdata('./NUVU/bin/dark.fits')
            dark_data = np.add(dark_data, newDark)
            print('i ' + str(newDark[7][7]))
        print('\n')
        dark_data = dark_data/100
        print('Overall: '+str(dark_data[7][7]))
        dark_data = np.delete(dark_data, 1, 0)
        dark_data = np.delete(dark_data, 0, 0)
        hdu = fits.PrimaryHDU(data=dark_data)
        hdu.writeto('./FITS/avgDark.fits', overwrite=True)

        image_data = np.zeros((128, 128))
        self.update_plot(rot_dial, sc1, sc2)

    def pull_image(self, rot_dial, sc1, sc2):
        print('Plotting Image')
        if dark_sub:
            new_image_data = np.subtract(image_data, dark_data)
        else:
            new_image_data = image_data

        centx, centy, disp, errs = figCalcs().get_centroids(new_image_data)
        dCx, dCy = np.meshgrid(dCenters, dCenters)
        gLines = np.linspace(-0.5, 127.5, image_mode + 1)
        rot_dial.setValue(180 + 60 * theta)
        self.plot_errors(sc2, errs)

        #sc1.axes.axes.imshow(new_image_data, vmin=lowerSlide, vmax=upperSlide, cmap='gray')
        if taking_dark:
            sc1.axes.axes.imshow(dark_data, cmap='gray')
        else:
            sc1.axes.axes.imshow(new_image_data, cmap='gray')
        sc1.axes.axes.set_xticks(gLines, minor=False)
        sc1.axes.axes.set_yticks(gLines, minor=False)
        sc1.axes.axes.set_xlim([-0.5, 127.5])
        sc1.axes.axes.set_ylim([127.5, -0.5])
        sc1.axes.axes.grid(color='limegreen', linestyle='-', linewidth=2, alpha=gla)
        sc1.axes.axes.plot(dCx, dCy, color='limegreen', ls='', marker='.', markersize=3, alpha=gca)
        sc1.axes.axes.plot(dCx + centx - centMu, dCy + centy - centMu, color='r', ls='', marker='.', markersize=3,
                           alpha=fca)
        sc1.axes.axes.set_xticklabels([])
        sc1.axes.axes.set_yticklabels([])
        #saveImage().save(sc1)

        if image_mode == 1:
            self.plotSquare(sc1)

        return sc1, sc2

    def plot_errors(self, sc2, errs):
        global zcoeffs
        zcoeffs = [li_coeff*x+(1-li_coeff)*y for x, y in zip(zcoeffs, errs)]

        if taking_dark:
            ymin = -1.0
            ymax = 1.0
            sc2.axes.axes.bar(z20_xs, np.zeros(20), color='k', edgecolor='limegreen', zorder=3)
        else:
            bound = np.max([np.abs(np.min(zcoeffs)), np.abs(np.max(zcoeffs))]) * 1.2
            if np.abs(bound) < 0.001:
                ymin = -1.0
                ymax = 1.0
            else:
                ymin = -bound
                ymax = bound
            sc2.axes.axes.bar(z20_xs, zcoeffs, color='k', edgecolor='limegreen', zorder=3)

        back = Rectangle((-0.5, ymin), 21, ymax-ymin, angle=0.0, color='k')
        sc2.axes.axes.add_patch(back)
        sc2.axes.axes.set_xticks(np.linspace(1, 20, 20), minor=False)
        sc2.axes.axes.set_yticks(np.linspace(ymin, ymax, 11), minor=False)
        sc2.axes.axes.grid(color='limegreen', linestyle='--', linewidth=0.5, which='both', zorder=0)
        sc2.axes.axes.axhline(y=0, color='limegreen')

        plt.xticks(fontsize=6, color='w')
        plt.yticks(fontsize=6, color='w')
        sc2.axes.axes.set_xlabel('Zernike Mode', color='w', fontsize=8)
        sc2.axes.axes.set_ylabel('Coefficient (nm)', color='w', fontsize=8)
        sc2.axes.axes.set_xlim([0.5, 20.5])
        sc2.axes.axes.set_ylim([ymin, ymax])
        return sc2

    def update_plot(self, rot_dial, sc1, sc2):
        sc1.axes.cla()
        sc2.axes.cla()
        self.pull_image(rot_dial, sc1, sc2)
        sc1.draw()
        sc2.draw()

    def plotSquare(self, sc1):
        sqMin = 63.5 - sl/2
        sqMax = 63.5 + sl/2
        normMin = (sqMin + 0.5)/128
        normMax = (sqMax + 0.5)/128
        sc1.axes.axes.axvline(x=sqMin, ymin=normMin, ymax=normMax, color='limegreen', ls='-', lw=1)
        sc1.axes.axes.axvline(x=sqMax, ymin=normMin, ymax=normMax, color='limegreen', ls='-', lw=1)
        sc1.axes.axes.axhline(y=sqMin, xmin=normMin, xmax=normMax, color='limegreen', ls='-', lw=1)
        sc1.axes.axes.axhline(y=sqMax, xmin=normMin, xmax=normMax, color='limegreen', ls='-', lw=1)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    MainWindow()
    app.exec_()
    #app.aboutToQuit.connect(MainWindow.closeEvent)