function hash = hash_mat_file(data)
%HASH_MAT_FILE  CRC-32 fingerprint of loaded .mat data (pure MATLAB, no Java).
%   hash = hash_mat_file(S) serializes S (e.g. the struct returned by load)
%   to bytes and returns its CRC-32 as an 8-character lowercase hex string.
%   Uses the standard IEEE 802.3 CRC-32 (same as zlib, PNG, zip).
%
%   A uint8 input is hashed as raw bytes, which allows a quick self-check:
%       hash_mat_file(uint8('123456789'))   % must return 'cbf43926'
%
%   Note: CRC-32 detects accidental changes to the data; it is not a
%   cryptographic hash.

    if isa(data, 'uint8')
        bytes = data(:);
    else
        bytes = getByteStreamFromArray(data);   % serialize struct to uint8
        bytes = bytes(:);
    end

    crc  = crc32_bytes(bytes);
    hash = lower(dec2hex(double(crc), 8));
end

function crc = crc32_bytes(bytes)
% Slicing-by-8 CRC-32: processes 8 bytes per loop iteration.

    persistent T
    m = uint32(255);

    if isempty(T)
        % T(:,1) is the standard byte table; T(:,k) extends it for slicing.
        T = zeros(256, 8, 'uint32');
        poly = uint32(hex2dec('EDB88320'));     % reflected CRC-32 polynomial
        for i = 0:255
            c = uint32(i);
            for j = 1:8
                if bitand(c, uint32(1))
                    c = bitxor(bitshift(c, -1), poly);
                else
                    c = bitshift(c, -1);
                end
            end
            T(i+1, 1) = c;
        end
        for k = 2:8
            prev = T(:, k-1);
            T(:, k) = bitxor(bitshift(prev, -8), T(double(bitand(prev, m)) + 1, 1));
        end
    end

    crc = uint32(4294967295);                   % 0xFFFFFFFF
    n   = numel(bytes);
    n8  = 8 * floor(n / 8);

    if n8 > 0
        w   = typecast(bytes(1:n8), 'uint32');  % little-endian 32-bit words
        one = w(1:2:end);
        two = w(2:2:end);

        % The second word of each 8-byte block doesn't depend on crc,
        % so its four table lookups are done for all blocks at once.
        X = bitxor( ...
                bitxor(T(double(bitand(two, m)) + 1, 4), ...
                       T(double(bitand(bitshift(two, -8), m)) + 1, 3)), ...
                bitxor(T(double(bitand(bitshift(two, -16), m)) + 1, 2), ...
                       T(double(bitshift(two, -24)) + 1, 1)));

        T5 = T(:,5); T6 = T(:,6); T7 = T(:,7); T8 = T(:,8);
        for k = 1:numel(one)
            a   = bitxor(crc, one(k));
            crc = bitxor( ...
                    bitxor(T8(double(bitand(a, m)) + 1), ...
                           T7(double(bitand(bitshift(a, -8), m)) + 1)), ...
                    bitxor(bitxor(T6(double(bitand(bitshift(a, -16), m)) + 1), ...
                                  T5(double(bitshift(a, -24)) + 1)), ...
                           X(k)));
        end
    end

    % Leftover bytes (fewer than 8), one at a time
    T1 = T(:,1);
    for k = n8+1:n
        crc = bitxor(bitshift(crc, -8), ...
                     T1(double(bitand(bitxor(crc, uint32(bytes(k))), m)) + 1));
    end

    crc = bitxor(crc, uint32(4294967295));
end
