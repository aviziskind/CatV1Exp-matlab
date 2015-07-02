function tmp_printDoubleRecords


    Gids_double = [227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 271 272 273 274 275 276 277 278 279 280 281 282 283 284 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299 300 301 302 303 304 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323 324 325 326 327 328 329 330 331 332 333 334 335 336 337 338 339 340 341 342 2196 2198 2200 2202 2204 2206 2214 2216 2218 2220 2224 2226 2228 2230 2234 2236 2238 2240 2242 2244 2246 2248 2250 2252 2254 2256 2258 2260 2262 2264 2268 2270 2272 2274 2276 2278 2282 2284 2286 2288 2290 2292 2294 2296 2298 2300 2302 2304 2306 2308 2310 2312 2314 2316 2318 2320 2322 2324 2326 2328 2330 2332 2334 2336 2338 2340 2342 2344 2346 2348 2350 2352 2354 2356 2358 2360 2362 2368 2370 2372 2374 2376 2378 2396 2398 2414 2416 2418 2420 2422 2424 2426 2428 2430 2432 2434 2436 2438 2440 2442 2444 2446 2448 2450 2452 2454 2456 2458 2460 2462 2464 2466 2468 2472 2474 2476 2478 2480 2482 2484 2486 2488 2490 2492 2494 2496 2498 2500 2502 2504 2506 2508 2510 2512 2514 2516 2518 2520 2522 2524 2526 2530 2532 2536 2538 2540 2542 2544 2546 2548 2550 4187 4189 4191 4224 4226 4228 4230 4232 4234 4236 4238 4240 4242 4244 4246 5306 5308];
    N = length(Gids_double);
    time_sec_diff_p = zeros(1, N); 
    time_sec_diff_l = zeros(1, N); 
    progressBar('init-', N)
    for gi = 1:N
%         progressBar(gi);
        Gid = Gids_double(gi);
        [hemi_L, hemi_R, time_sec_diff_p(gi), time_sec_diff_l(gi)] = printGidRecords(Gid);
        fprintf('%s (pen) %s (loc)\n', sec2hms(time_sec_diff_p(gi)), sec2hms(time_sec_diff_l(gi)));
        
        assert(length(unique(hemi_L))==1);
        assert(length(unique(hemi_R))==1);        
        
        assert(all( xor(hemi_L, hemi_R)) )
        3;
                
    end
    3;

end
            
function [hemi_L, hemi_R, diff_time_pen_sec, diff_time_loc_sec] = printGidRecords(Gid)


    T1 = {'TBL_AP_ML_ZERO', 'ELECTRODE_ID', 'TBL_ELECTRODES'};
    T2 = {T1, 'TBL_AP_ML_ZERO', 'AP_ML_ZERO_ID', 'TBL_PENETRATIONS'};
    T3 = {'TBL_ANIMALS', 'ANIMAL_ID', T2, 'TBL_ELECTRODES'};
    T4 = {T3, 'TBL_PENETRATIONS', 'PENETRATION_ID', 'TBL_LOCATIONS'};
    T5 = {'TBL_DATA_FILES', 'DATAFILE_ID', 'TBL_LOCS_FILES_LINKS'};
    T6 = {'TBL_LOC_DEPTHS', 'LOC_DEPTH_ID',  T5, 'TBL_LOCS_FILES_LINKS'};
    T7 = {T4, 'TBL_LOCATIONS', 'LOCATION_ID', T6, 'TBL_LOC_DEPTHS'};
    T8 = {T7, 'TBL_DATA_FILES', 'DATAFILE_ID', 'TBL_GROUPS'};
    T9 = {'TBL_ELECTRODE_TYPES', 'ELECTRODE_TYPE_ID', T8, 'TBL_ELECTRODES'};
    joinedTables = T9;

    

    hnd = dbOpenExpDb;
    criterea = {'TBL_GROUPS.GROUP_ID', Gid};
    [recs, fld_names] = getFieldsFromDatabaseTable(hnd, '*', joinedTables, criterea);
    rec_diff = zeros(size(fld_names));
    for i = 1:length(recs)
        rec_diff(i) = ~isequalwithequalnans(recs(1,i), recs(2,i));
    end
    [fld_names, recs', num2cell(rec_diff)];
    
    hemi_l_idx = find(strcmp(fld_names, 'BLN_HEMISPHERE_L'), 1);
    hemi_r_idx = find(strcmp(fld_names, 'BLN_HEMISPHERE_R'), 1);
    
    hemi_L = [recs{:,hemi_l_idx}];
    hemi_R = [recs{:,hemi_r_idx}];
    3;
    
    
    pen_time_idx = find(strcmp(fld_names, 'MEM_PENETRATION_NOTES'), 1)+1;  % DTM_CREATED is just after this one    
    loc_time_idx = find(strcmp(fld_names, 'BLN_NO_DATA_FILES'), 1)+1;  % DTM_CREATED is just after this one
    
    pen_time_str = [recs(:,pen_time_idx)];
    loc_time_str = [recs(:,loc_time_idx)];
    
    dateMask = 'mm/dd/yyyy HH:MM:SS PM';
    date2num = @(date_str) datenum(date_str, dateMask);
    pen_time = cellfun(date2num, pen_time_str);
    loc_time = cellfun(date2num, loc_time_str);

    secPerDay = 24*60*60;
    diff_time_pen_sec = diff(pen_time)*secPerDay;
    diff_time_loc_sec = diff(loc_time)*secPerDay;
    3;
    
    nTrunc = 40;
%     clc;
return;
    for i = 1:length(recs)
        if strcmp(fld_names{i}, 'MEM_SPIKETIMES_MTX')
            3;
        end
        rec_diff_str = iff(rec_diff(i), '<=******', '');
        
        fprintf('%20s : %40s, %40s, %s\n', trunc(fld_names{i}, nTrunc), trunc(var2str(recs{1,i}), nTrunc), trunc(var2str(recs{2,i}), nTrunc), rec_diff_str)                 
    end
3;

end            

function s = trunc(s, n)
    if length(s) > n
        s = s(1:n);
    end        
    s(s==13) = '/';
    s(s==10) = '';
end
