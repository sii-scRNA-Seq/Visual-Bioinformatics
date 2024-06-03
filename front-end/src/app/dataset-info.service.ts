import { BehaviorSubject, Observable } from 'rxjs';
import { Injectable } from '@angular/core';

import { BackendHttpClient } from './backend-http.client';
import { DatasetInfoServiceInterface } from './dataset-info.service.interface';

@Injectable({
  providedIn: 'root'
})
export class DatasetInfoService implements DatasetInfoServiceInterface {
  private readonly datasetInfo$: BehaviorSubject<string[]> = new BehaviorSubject<string[]> ([]);
  readonly datasetInfo: Observable<string[]> = this.datasetInfo$.asObservable();

  constructor(private backendHttpClient: BackendHttpClient) { 
    this.setDatasetInfo();
  }

  async setDatasetInfo(): Promise<void> {
    const datasetInfo = await this.backendHttpClient.getDatasetInfo();
    this.datasetInfo$.next(datasetInfo);
  }
}
