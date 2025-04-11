function segments = Fun_SegmentContour(corner_idx, n)
%% generate a cell-type "segments" which are several segments splitted by corner_idx in region [1:n]
    % corner_idx: the splitted location
    % n: determine the range of the interval
    % segments: the splitted segments
%%
    segments = {};
    prev = corner_idx(1);
    for i = 2:length(corner_idx)
        next = corner_idx(i);
        segments{end+1} = prev:next;
        prev = next;
    end
    segments{end+1} = [prev:n,1:corner_idx(1)-1];
end
