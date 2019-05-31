function m = colorinput(pic)
% This function returns the COLOR DOUBLE OPPONENT CELL responses
% for an rgb-image. 
% The argument pic should be a rgb-color-image
% (size(pic) == [x y 3]), 
% It returns a cell-array (4x size(pic,1),size(pic,2))
% where m{1} is red-center-green-surround
%       m{2} is green-center-red-surround
%       m{3} is blue-center-yellow-surround
%   and m{4} is yellow-center-blue-surround cell input.
%
% example:
%    pic = imread('horse.jpg');
%    m = colorinput(pic);
%    for i=1:4,
%      subplot(2,2,i), imagesc(m{i}); colorbar;
%    end; colormap(gray);
% (Hauke)

%produce a dog filter
fs = 31; %% filter size in pixels
f = zeros(fs, fs);
sigmae = 3.6; sigmai = 3.8;
for i=1:fs,
  for j=1:fs,
    x = i-round(fs/2)+1;
    y = j-round(fs/2)+1;
    f(i,j) = 1/(sqrt(2*pi)*sigmae) * exp(-(x^2 + y^2)/(2*sigmae^2))-...
        1/(sqrt(2*pi)*sigmai) * exp(-(x^2+y^2)/(2*sigmai^2));
  end;
end;
surfl(f); colormap(gray);

% to be sure that we have doubles not uint8
pic = double(pic); 

redch(:,:) = pic(:,:,1);
greench(:,:) = pic(:,:,2);
bluech(:,:) = pic(:,:,3);

% red-center-green-surround opponent cell output
RGC = filter2(f, redch,'valid') + filter2(-f, greench, 'valid');

% green-center-red-surround opponent cell
GRC = filter2(f, greench, 'valid') + filter2(-f, redch, 'valid');

% blue-center-yellow-surround opponent cell
BYC = filter2(f, bluech, 'valid') + filter2(-f, redch+greench, 'valid');

% yellow-center-blue-surround opponent cell
YBC = filter2(f, (redch+greench), 'valid') + filter2(-f, bluech, 'valid');

figure, colormap(gray);
subplot(2,2,1), imagesc(RGC); title('RGC'); colorbar;
subplot(2,2,2), imagesc(GRC); title('GRC');colorbar;
subplot(2,2,3), imagesc(BYC); title('BYC');colorbar;
subplot(2,2,4), imagesc(YBC); title('YBC');colorbar;
colormap(gray);

m{1} = RGC; m{2} = GRC; m{3} = BYC; m{4} = YBC;

