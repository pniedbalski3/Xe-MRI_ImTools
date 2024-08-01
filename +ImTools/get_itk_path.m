function itk_folder = get_itk_path()

sys_folds = dir('C:\Program Files');

sys_folds = struct2cell(sys_folds);

mynames = sys_folds(1,:);

itk = find(contains(mynames,'ITK-SNAP'));

itk_folder = mynames{itk};

