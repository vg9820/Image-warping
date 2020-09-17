%% read image
clear all
im = imread('C:\Users\wzw98\Desktop\���ѵ��\Matlab-1_Assignment\TestImage\warp_test.bmp');

%% draw 2 copies of the image
figure('Units', 'pixel', 'Position', [100,100,1000,700], 'toolbar', 'none');
subplot(121); imshow(im); title({'Source image', 'Press the red tool button to add point-point constraints'});
subplot(122); himg = imshow(im*0); title({'Warpped Image', 'Press the blue tool button to compute the warpped image'});
hToolPoint = uipushtool('CData', reshape(repmat([1 0 0], 100, 1), [10 10 3]), 'TooltipString', 'add point constraints to the map', ...
                        'ClickedCallback', @toolPositionCB, 'UserData', []);
hToolWarp = uipushtool('CData', reshape(repmat([0 0 1], 100, 1), [10 10 3]), 'TooltipString', 'compute warped image', ...
                       'ClickedCallback', @toolWarpCB);

%% TODO: implement function: RBFImageWarp
% check the title above the image for how to use the simple user-interface to define point-constraints and compute the warpped image