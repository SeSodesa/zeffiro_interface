function zef = zef_dataBank_refreshTree(zef)

if nargin == 0
zef = evalin('base','zef')
end

zef.dataBank.app.Tree.Children.delete;
zef_dataBank_hash2tree;

if nargout == 0
assignin('base','zef',zef);
end

end
