import { first, firstValueFrom } from 'rxjs';
import { HttpClientTestingModule } from '@angular/common/http/testing';
import { TestBed } from '@angular/core/testing';

import { BackendHttpClient } from './backend-http.client';
import { DatasetInfoService } from './dataset-info.service';

describe('DatasetInfoService', () => {
  let service: DatasetInfoService;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [
        HttpClientTestingModule,
      ]
    });
    service = TestBed.inject(DatasetInfoService);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });

  describe('setDatasetInfo', () => {
    it('should call client getDatasetInfo method', async () => {
      const backendHttpClient: BackendHttpClient = TestBed.inject(BackendHttpClient);
      const backendHttpClientSpy = spyOn(backendHttpClient, 'getDatasetInfo');
      await service.setDatasetInfo();
      expect(backendHttpClientSpy).toHaveBeenCalled();
    });

    it('should update datasetInfo attribute as expected', async () => {
      const backendHttpClient: BackendHttpClient = TestBed.inject(BackendHttpClient);
      spyOn(backendHttpClient, 'getDatasetInfo').and.returnValue(Promise.resolve([
        {key: 'option1', title: 'Option 1'},
        {key: 'option2', title: 'Option 2'}
      ]));
      await service.setDatasetInfo();
      const response = await firstValueFrom(service.datasetInfo.pipe(first()));
      expect(response).toEqual([
        {key: 'option1', title: 'Option 1'},
        {key: 'option2', title: 'Option 2'}
      ]);
    });
  });
});
