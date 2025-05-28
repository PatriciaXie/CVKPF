classdef ImageReader
    properties
        nCamera
        datasetPath
        fileType
        filename
        startFrame
        endFrame
    end
    
    methods
        function imgHandle = ImageReader(nCamera, datasetPath, fileType)
            % ImageReader ��ȡ����ͼͼ������
            imgHandle.nCamera = nCamera;
            imgHandle.datasetPath = datasetPath;
            imgHandle.fileType = fileType;
            for v = 1:1:nCamera
                info = dir(fullfile(sprintf('%s/cam%d', datasetPath, v), sprintf('*.%s', fileType)));
                if isempty(info)
                    error('%s�ļ���û��ͼ������', sprintf('%s/cam%d', datasetPath, v));
                else
                    imgHandle.filename{v} = {info.name}';
                end
            end
            imgHandle.startFrame = 1;
            imgHandle.endFrame = length(imgHandle.filename{1});
        end
        
        function img = Read(imgHandle, idxCamera, idxFrames)
            try
                idxFrame = idxFrames(1);
                file =  sprintf('%s/cam%d/%s', imgHandle.datasetPath, idxCamera, imgHandle.filename{idxCamera}{idxFrame});
                img = imread(file);
                img = double(img)/255;
                if ~isscalar(idxFrames)
                    nImg = length(idxFrames);
                    imgs = zeros([size(img), nImg]);
                    imgs(:, :, 1) = img;
                    for i = 2:1:nImg
                        idxFrame = idxFrames(i);
                        file =  sprintf('%s/cam%d/%s', imgHandle.datasetPath, idxCamera, imgHandle.filename{idxCamera}{idxFrame});
                        img = imread(file);
                        img = double(img)/255;
                        imgs(:, :, i) = img;
                    end
                    img = imgs;
                end
            catch
                error('ͼƬ��ȡ���󳬳����ݿⷶΧ������Camera%d-Frame%d�����ݷ�ΧFrame%d-%d', idxCamera, idxFrame, imgHandle.startFrame, imgHandle.endFrame);
            end
        end
    end
            
end