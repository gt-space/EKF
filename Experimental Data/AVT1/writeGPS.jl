using AstroTime


using CSV
using Tables


gpstxt = open(normpath(joinpath((@__FILE__), raw"..\gps_overnight.txt")), "r")
gpslines = readlines(gpstxt)


initial_index = 10
final_index = 33806
i = 1
j = initial_index
φ = zeros(8450)
λ = zeros(8450)
h = zeros(8450)
vx = zeros(8450)
vy = zeros(8450)
vz = zeros(8450)
t = zeros(8450)


while true
    temp_string = gpslines[j]
    temp_string = chop(temp_string, head = 14, tail = 1)
    gps_state = split(temp_string, ['°'])
    φ[i] = parse(Float64, gps_state[1])
    λ[i] = parse(Float64, chop(gps_state[2], head = 6, tail = 0))
    h[i] = parse(Float64, chop(gps_state[3], head = 6, tail = 0))
    temp_string = gpslines[j + 1]
    temp_string = chop(temp_string, head = 22, tail = 4)
    velocity_state = split(temp_string, [','])
    vx[i] = parse(Float64, chop(velocity_state[1], head = 0, tail = 4))
    vy[i] = parse(Float64, chop(velocity_state[2], head = 6, tail = 4))
    vz[i] = parse(Float64, chop(velocity_state[3], head = 6, tail = 0))
    temp_string = gpslines[j + 2]
    temp_string = chop(temp_string, head = 6, tail = 4)
    first_split = split(temp_string, [' '])
    date = split(first_split[1], ['-'])
    time = split(first_split[2], [':'])
    t[i] = value(seconds(Epoch{UniversalTime}(parse(Int, date[1]), parse(Int, date[2]), parse(Int, date[3]), parse(Int, time[1]), parse(Int, time[2]), parse(Float64, time[3])) - Epoch{UniversalTime}(2000, 1, 1, 0, 0, 0.0)))
    if j == final_index
        break
    end
    i += 1
    j += 4
end
t .-= t[1]


my_mat = [φ λ h vx vy vz t]
CSV.write(normpath(joinpath((@__FILE__), raw"..\GPSdata.csv")), Tables.table(my_mat), header = ["Geodetic Latitude", "Geodetic Longitude", "Geodetic Altitude", "x-velocity", "y-velocity", "z-velocity", "time from start"])







