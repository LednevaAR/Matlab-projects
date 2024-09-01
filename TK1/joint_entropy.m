function H = joint_entropy(X, Y)
    H = entropy(string(strcat(X.', Y.')).');
end

