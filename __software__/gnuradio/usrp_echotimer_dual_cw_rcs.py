#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: Usrp Echotimer Dual Cw Rcs
# GNU Radio version: 3.8.2.0

from distutils.version import StrictVersion

if __name__ == '__main__':
    import ctypes
    import sys
    if sys.platform.startswith('linux'):
        try:
            x11 = ctypes.cdll.LoadLibrary('libX11.so')
            x11.XInitThreads()
        except:
            print("Warning: failed to XInitThreads()")

from gnuradio import blocks
from gnuradio import filter
from gnuradio.filter import firdes
from gnuradio import gr
import sys
import signal
from PyQt5 import Qt
from argparse import ArgumentParser
from gnuradio.eng_arg import eng_float, intx
from gnuradio import eng_notation
from gnuradio.qtgui import Range, RangeWidget
import radar

from gnuradio import qtgui

class usrp_echotimer_dual_cw_rcs(gr.top_block, Qt.QWidget):

    def __init__(self):
        gr.top_block.__init__(self, "Usrp Echotimer Dual Cw Rcs")
        Qt.QWidget.__init__(self)
        self.setWindowTitle("Usrp Echotimer Dual Cw Rcs")
        qtgui.util.check_set_qss()
        try:
            self.setWindowIcon(Qt.QIcon.fromTheme('gnuradio-grc'))
        except:
            pass
        self.top_scroll_layout = Qt.QVBoxLayout()
        self.setLayout(self.top_scroll_layout)
        self.top_scroll = Qt.QScrollArea()
        self.top_scroll.setFrameStyle(Qt.QFrame.NoFrame)
        self.top_scroll_layout.addWidget(self.top_scroll)
        self.top_scroll.setWidgetResizable(True)
        self.top_widget = Qt.QWidget()
        self.top_scroll.setWidget(self.top_widget)
        self.top_layout = Qt.QVBoxLayout(self.top_widget)
        self.top_grid_layout = Qt.QGridLayout()
        self.top_layout.addLayout(self.top_grid_layout)

        self.settings = Qt.QSettings("GNU Radio", "usrp_echotimer_dual_cw_rcs")

        try:
            if StrictVersion(Qt.qVersion()) < StrictVersion("5.0.0"):
                self.restoreGeometry(self.settings.value("geometry").toByteArray())
            else:
                self.restoreGeometry(self.settings.value("geometry"))
        except:
            pass

        ##################################################
        # Variables
        ##################################################
        self.samp_rate = samp_rate = 14250000
        self.packet_len = packet_len = 2**21
        self.freq_res = freq_res = samp_rate/float(packet_len)
        self.freq = freq = (-3000000,3000000)
        self.center_freq = center_freq = 2.4e9
        self.v_res = v_res = freq_res*3e8/2/center_freq
        self.time_res = time_res = packet_len/float(samp_rate)
        self.threshold = threshold = -200
        self.samp_protect = samp_protect = 1
        self.range_time = range_time = 60
        self.range_res = range_res = (3e8/2/float((freq[1]-freq[0])))
        self.min_output_buffer = min_output_buffer = int(packet_len*2)
        self.max_output_buffer = max_output_buffer = 0
        self.gain_tx = gain_tx = 40
        self.gain_rx = gain_rx = 10
        self.delay_samp = delay_samp = 33
        self.decim_fac = decim_fac = 2**7
        self.amplitude = amplitude = 0.5

        ##################################################
        # Blocks
        ##################################################
        self._threshold_range = Range(-500, 100, 1, -200, 200)
        self._threshold_win = RangeWidget(self._threshold_range, self.set_threshold, 'Find peak threshold', "counter_slider", float)
        self.top_grid_layout.addWidget(self._threshold_win, 1, 0, 1, 1)
        for r in range(1, 2):
            self.top_grid_layout.setRowStretch(r, 1)
        for c in range(0, 1):
            self.top_grid_layout.setColumnStretch(c, 1)
        self._samp_protect_range = Range(0, 100, 1, 1, 200)
        self._samp_protect_win = RangeWidget(self._samp_protect_range, self.set_samp_protect, 'Find peak protected samples', "counter_slider", float)
        self.top_grid_layout.addWidget(self._samp_protect_win, 1, 1, 1, 1)
        for r in range(1, 2):
            self.top_grid_layout.setRowStretch(r, 1)
        for c in range(1, 2):
            self.top_grid_layout.setColumnStretch(c, 1)
        self._gain_tx_range = Range(0, 100, 1, 40, 200)
        self._gain_tx_win = RangeWidget(self._gain_tx_range, self.set_gain_tx, 'TX gain', "counter_slider", float)
        self.top_grid_layout.addWidget(self._gain_tx_win, 0, 0, 1, 1)
        for r in range(0, 1):
            self.top_grid_layout.setRowStretch(r, 1)
        for c in range(0, 1):
            self.top_grid_layout.setColumnStretch(c, 1)
        self._gain_rx_range = Range(0, 100, 1, 10, 200)
        self._gain_rx_win = RangeWidget(self._gain_rx_range, self.set_gain_rx, 'RX gain', "counter_slider", float)
        self.top_grid_layout.addWidget(self._gain_rx_win, 0, 1, 1, 1)
        for r in range(0, 1):
            self.top_grid_layout.setRowStretch(r, 1)
        for c in range(1, 2):
            self.top_grid_layout.setColumnStretch(c, 1)
        self._delay_samp_range = Range(0, 100, 1, 33, 200)
        self._delay_samp_win = RangeWidget(self._delay_samp_range, self.set_delay_samp, 'Number delay samples', "counter_slider", float)
        self.top_grid_layout.addWidget(self._delay_samp_win)
        self.rational_resampler_xxx_0_0 = filter.rational_resampler_ccc(
                interpolation=1,
                decimation=decim_fac,
                taps=None,
                fractional_bw=None)
        self.rational_resampler_xxx_0 = filter.rational_resampler_ccc(
                interpolation=1,
                decimation=decim_fac,
                taps=None,
                fractional_bw=None)
        self.radar_usrp_echotimer_cc_0 = radar.usrp_echotimer_cc(samp_rate, center_freq, int(delay_samp), 'addr=10.10.3.20', '', 'internal', 'none', 'TX/RX', gain_tx, 0.1, 0.05, 0, 'addr=10.10.3.25', '', 'mimo', 'mimo', 'TX/RX', gain_rx, 0.1, 0.05, 0, "packet_len")
        self.radar_ts_fft_cc_0_0 = radar.ts_fft_cc(packet_len//decim_fac,  "packet_len")
        self.radar_ts_fft_cc_0_0.set_min_output_buffer(4194304)
        self.radar_ts_fft_cc_0 = radar.ts_fft_cc(packet_len//decim_fac,  "packet_len")
        self.radar_ts_fft_cc_0.set_min_output_buffer(4194304)
        self.radar_signal_generator_cw_c_0_0 = radar.signal_generator_cw_c(packet_len, samp_rate, [freq[1]], amplitude, "packet_len")
        self.radar_signal_generator_cw_c_0_0.set_min_output_buffer(4194304)
        self.radar_signal_generator_cw_c_0 = radar.signal_generator_cw_c(packet_len, samp_rate, [freq[0]], amplitude, "packet_len")
        self.radar_signal_generator_cw_c_0.set_min_output_buffer(4194304)
        self.radar_qtgui_time_plot_0 = radar.qtgui_time_plot(100, 'range', [0,range_res], range_time, '')
        self.radar_find_max_peak_c_0 = radar.find_max_peak_c(samp_rate//decim_fac, threshold, int(samp_protect), [-300,300], True, "packet_len")
        self.radar_estimator_fsk_0 = radar.estimator_fsk(center_freq, (freq[1]-freq[0]), False)
        self.blocks_tagged_stream_multiply_length_0_0 = blocks.tagged_stream_multiply_length(gr.sizeof_gr_complex*1, "packet_len", 1.0/float(decim_fac))
        self.blocks_tagged_stream_multiply_length_0_0.set_min_output_buffer(4194304)
        self.blocks_tagged_stream_multiply_length_0 = blocks.tagged_stream_multiply_length(gr.sizeof_gr_complex*1, "packet_len", 1.0/float(decim_fac))
        self.blocks_tagged_stream_multiply_length_0.set_min_output_buffer(4194304)
        self.blocks_multiply_conjugate_cc_1 = blocks.multiply_conjugate_cc(1)
        self.blocks_multiply_conjugate_cc_1.set_min_output_buffer(4194304)
        self.blocks_multiply_conjugate_cc_0_0 = blocks.multiply_conjugate_cc(1)
        self.blocks_multiply_conjugate_cc_0_0.set_min_output_buffer(4194304)
        self.blocks_multiply_conjugate_cc_0 = blocks.multiply_conjugate_cc(1)
        self.blocks_multiply_conjugate_cc_0.set_min_output_buffer(4194304)
        self.blocks_add_xx_1 = blocks.add_vcc(1)
        self.blocks_add_xx_1.set_min_output_buffer(4194304)



        ##################################################
        # Connections
        ##################################################
        self.msg_connect((self.radar_estimator_fsk_0, 'Msg out'), (self.radar_qtgui_time_plot_0, 'Msg in'))
        self.msg_connect((self.radar_find_max_peak_c_0, 'Msg out'), (self.radar_estimator_fsk_0, 'Msg in'))
        self.connect((self.blocks_add_xx_1, 0), (self.radar_usrp_echotimer_cc_0, 0))
        self.connect((self.blocks_multiply_conjugate_cc_0, 0), (self.rational_resampler_xxx_0, 0))
        self.connect((self.blocks_multiply_conjugate_cc_0_0, 0), (self.rational_resampler_xxx_0_0, 0))
        self.connect((self.blocks_multiply_conjugate_cc_1, 0), (self.radar_find_max_peak_c_0, 0))
        self.connect((self.blocks_tagged_stream_multiply_length_0, 0), (self.radar_ts_fft_cc_0, 0))
        self.connect((self.blocks_tagged_stream_multiply_length_0_0, 0), (self.radar_ts_fft_cc_0_0, 0))
        self.connect((self.radar_signal_generator_cw_c_0, 0), (self.blocks_add_xx_1, 0))
        self.connect((self.radar_signal_generator_cw_c_0, 0), (self.blocks_multiply_conjugate_cc_0, 1))
        self.connect((self.radar_signal_generator_cw_c_0_0, 0), (self.blocks_add_xx_1, 1))
        self.connect((self.radar_signal_generator_cw_c_0_0, 0), (self.blocks_multiply_conjugate_cc_0_0, 1))
        self.connect((self.radar_ts_fft_cc_0, 0), (self.blocks_multiply_conjugate_cc_1, 0))
        self.connect((self.radar_ts_fft_cc_0_0, 0), (self.blocks_multiply_conjugate_cc_1, 1))
        self.connect((self.radar_usrp_echotimer_cc_0, 0), (self.blocks_multiply_conjugate_cc_0, 0))
        self.connect((self.radar_usrp_echotimer_cc_0, 0), (self.blocks_multiply_conjugate_cc_0_0, 0))
        self.connect((self.rational_resampler_xxx_0, 0), (self.blocks_tagged_stream_multiply_length_0, 0))
        self.connect((self.rational_resampler_xxx_0_0, 0), (self.blocks_tagged_stream_multiply_length_0_0, 0))


    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "usrp_echotimer_dual_cw_rcs")
        self.settings.setValue("geometry", self.saveGeometry())
        event.accept()

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.set_freq_res(self.samp_rate/float(self.packet_len))
        self.set_time_res(self.packet_len/float(self.samp_rate))

    def get_packet_len(self):
        return self.packet_len

    def set_packet_len(self, packet_len):
        self.packet_len = packet_len
        self.set_freq_res(self.samp_rate/float(self.packet_len))
        self.set_min_output_buffer(int(self.packet_len*2))
        self.set_time_res(self.packet_len/float(self.samp_rate))

    def get_freq_res(self):
        return self.freq_res

    def set_freq_res(self, freq_res):
        self.freq_res = freq_res
        self.set_v_res(self.freq_res*3e8/2/self.center_freq)

    def get_freq(self):
        return self.freq

    def set_freq(self, freq):
        self.freq = freq
        self.set_range_res((3e8/2/float((self.freq[1]-self.freq[0]))))

    def get_center_freq(self):
        return self.center_freq

    def set_center_freq(self, center_freq):
        self.center_freq = center_freq
        self.set_v_res(self.freq_res*3e8/2/self.center_freq)

    def get_v_res(self):
        return self.v_res

    def set_v_res(self, v_res):
        self.v_res = v_res

    def get_time_res(self):
        return self.time_res

    def set_time_res(self, time_res):
        self.time_res = time_res

    def get_threshold(self):
        return self.threshold

    def set_threshold(self, threshold):
        self.threshold = threshold
        self.radar_find_max_peak_c_0.set_threshold(self.threshold)

    def get_samp_protect(self):
        return self.samp_protect

    def set_samp_protect(self, samp_protect):
        self.samp_protect = samp_protect
        self.radar_find_max_peak_c_0.set_samp_protect(int(self.samp_protect))

    def get_range_time(self):
        return self.range_time

    def set_range_time(self, range_time):
        self.range_time = range_time

    def get_range_res(self):
        return self.range_res

    def set_range_res(self, range_res):
        self.range_res = range_res

    def get_min_output_buffer(self):
        return self.min_output_buffer

    def set_min_output_buffer(self, min_output_buffer):
        self.min_output_buffer = min_output_buffer

    def get_max_output_buffer(self):
        return self.max_output_buffer

    def set_max_output_buffer(self, max_output_buffer):
        self.max_output_buffer = max_output_buffer

    def get_gain_tx(self):
        return self.gain_tx

    def set_gain_tx(self, gain_tx):
        self.gain_tx = gain_tx
        self.radar_usrp_echotimer_cc_0.set_tx_gain(self.gain_tx)

    def get_gain_rx(self):
        return self.gain_rx

    def set_gain_rx(self, gain_rx):
        self.gain_rx = gain_rx
        self.radar_usrp_echotimer_cc_0.set_rx_gain(self.gain_rx)

    def get_delay_samp(self):
        return self.delay_samp

    def set_delay_samp(self, delay_samp):
        self.delay_samp = delay_samp
        self.radar_usrp_echotimer_cc_0.set_num_delay_samps(int(self.delay_samp))

    def get_decim_fac(self):
        return self.decim_fac

    def set_decim_fac(self, decim_fac):
        self.decim_fac = decim_fac
        self.blocks_tagged_stream_multiply_length_0.set_scalar(1.0/float(self.decim_fac))
        self.blocks_tagged_stream_multiply_length_0_0.set_scalar(1.0/float(self.decim_fac))

    def get_amplitude(self):
        return self.amplitude

    def set_amplitude(self, amplitude):
        self.amplitude = amplitude





def main(top_block_cls=usrp_echotimer_dual_cw_rcs, options=None):

    if StrictVersion("4.5.0") <= StrictVersion(Qt.qVersion()) < StrictVersion("5.0.0"):
        style = gr.prefs().get_string('qtgui', 'style', 'raster')
        Qt.QApplication.setGraphicsSystem(style)
    qapp = Qt.QApplication(sys.argv)

    tb = top_block_cls()

    tb.start()

    tb.show()

    def sig_handler(sig=None, frame=None):
        Qt.QApplication.quit()

    signal.signal(signal.SIGINT, sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)

    timer = Qt.QTimer()
    timer.start(500)
    timer.timeout.connect(lambda: None)

    def quitting():
        tb.stop()
        tb.wait()

    qapp.aboutToQuit.connect(quitting)
    qapp.exec_()

if __name__ == '__main__':
    main()
