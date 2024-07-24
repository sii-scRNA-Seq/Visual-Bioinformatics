import { first, firstValueFrom } from 'rxjs';
import { provideHttpClientTesting } from '@angular/common/http/testing';
import { TestBed } from '@angular/core/testing';

import { BackendHttpClient } from './backend-http.client';
import { DatasetInfoService } from './dataset-info.service';
import { provideHttpClient, withInterceptorsFromDi } from '@angular/common/http';

describe('DatasetInfoService', () => {
  let service: DatasetInfoService;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [],
      providers: [provideHttpClient(withInterceptorsFromDi()), provideHttpClientTesting()]
    });
  });

  it('should be created', () => {
    service = TestBed.inject(DatasetInfoService);
    expect(service).toBeTruthy();
  });

  describe('constructor', () => {
    it('should result in datasetInfo being loaded', async () => {
      const backendHttpClient: BackendHttpClient = TestBed.inject(BackendHttpClient);
      spyOn(backendHttpClient, 'getDatasetInfo').and.returnValue(Promise.resolve([
        {key: 'option1', title: 'Option 1', integration_obs: []},
        {key: 'option2', title: 'Option 2', integration_obs: []}
      ]));
      service = await TestBed.inject(DatasetInfoService);
      const datasetInfo = await firstValueFrom(service.datasetInfo.pipe(first()));
      expect(datasetInfo).toEqual([
        {key: 'option1', title: 'Option 1', integration_obs: []},
        {key: 'option2', title: 'Option 2', integration_obs: []}
      ]);
    });
  });

  describe('setDatasetInfo', () => {
    beforeEach(() => {
      service = TestBed.inject(DatasetInfoService);
    });

    it('should call client getDatasetInfo method', async () => {
      const backendHttpClient: BackendHttpClient = TestBed.inject(BackendHttpClient);
      const backendHttpClientSpy = spyOn(backendHttpClient, 'getDatasetInfo');
      await service.setDatasetInfo();
      expect(backendHttpClientSpy).toHaveBeenCalled();
    });

    it('should update datasetInfo attribute as expected', async () => {
      const backendHttpClient: BackendHttpClient = TestBed.inject(BackendHttpClient);
      spyOn(backendHttpClient, 'getDatasetInfo').and.returnValue(Promise.resolve([
        {key: 'option1', title: 'Option 1', integration_obs: []},
        {key: 'option2', title: 'Option 2', integration_obs: []}
      ]));
      await service.setDatasetInfo();
      const response = await firstValueFrom(service.datasetInfo.pipe(first()));
      expect(response).toEqual([
        {key: 'option1', title: 'Option 1', integration_obs: []},
        {key: 'option2', title: 'Option 2', integration_obs: []}
      ]);
    });
  });
});
