function out = tricks_common()
%% tricks_common - Shared utilities for Chapter 21 (Special Tricks).
%
% Returns a struct of helpers:
%   out.NAVY, .CORAL, .TEAL, .PURPLE, .GOLD, .OLIVE  (RGB triples in [0,1])
%   out.configure_style                              (function handle)
%   out.output_dir(script_dir)                       (function handle)
%
% Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
% Part of "Computational Etudes: A Spectral Approach"

    out.NAVY   = [0x14, 0x2D, 0x6E] / 255;
    out.SKY    = [0x78, 0x96, 0xD2] / 255;
    out.CORAL  = [0xE7, 0x4C, 0x3C] / 255;
    out.TEAL   = [0x16, 0xA0, 0x85] / 255;
    out.PURPLE = [0x8E, 0x44, 0xAD] / 255;
    out.ORANGE = [0xE6, 0x7E, 0x22] / 255;
    out.GOLD   = [0xD4, 0xA0, 0x17] / 255;
    out.OLIVE  = [0x6B, 0x8E, 0x23] / 255;
    out.configure_style = @configure_style;
    out.output_dir      = @output_dir;
end

function configure_style()
    set(groot, 'defaultAxesFontName', 'CMU Serif');
    set(groot, 'defaultTextFontName', 'CMU Serif');
    set(groot, 'defaultAxesFontSize', 11);
    set(groot, 'defaultAxesLineWidth', 0.8);
    set(groot, 'defaultLineLineWidth', 1.0);
    set(groot, 'defaultAxesBox', 'on');
end

function dir_out = output_dir(script_dir)
    p = script_dir;
    while ~isempty(p) && ~exist(fullfile(p, 'textbook'), 'dir')
        new_p = fileparts(p);
        if strcmp(new_p, p); break; end
        p = new_p;
    end
    dir_out = fullfile(p, 'textbook', 'figures', 'ch21', 'matlab');
    if ~exist(dir_out, 'dir'); mkdir(dir_out); end
end
