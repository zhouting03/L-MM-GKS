function lplq_write_csv(path,rows)
% 写标准 CSV，不依赖 table/writetable，兼容 MATLAB 与 Octave。
if isempty(rows), return; end
fid=fopen(path,'w'); if fid<0, error('LpLq:Output','Cannot create %s',path); end
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fields=fieldnames(rows{1}); fprintf(fid,'%s\n',strjoin(fields,','));
for j=1:numel(rows)
    values=cell(size(fields));
    for k=1:numel(fields)
        v=rows{j}.(fields{k});
        if isnumeric(v) || islogical(v), values{k}=sprintf('%.12g',v);
        else, values{k}=['"',strrep(char(v),'"','""'),'"']; end
    end
    fprintf(fid,'%s\n',strjoin(values,','));
end
end
