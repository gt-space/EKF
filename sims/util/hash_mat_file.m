function hash = hash_mat_file(filename)

    fid = fopen(filename, 'rb');

    if fid == -1
        error('Could not open file: %s', filename);
    end

    cleanup = onCleanup(@() fclose(fid));

    md = java.security.MessageDigest.getInstance('SHA-256');

    while true
        data = fread(fid, 1e6, '*uint8');

        if isempty(data)
            break;
        end

        md.update(data);
    end

    digest = typecast(md.digest(), 'uint8');
    hash = lower(reshape(dec2hex(digest, 2).', 1, []));

end